/* g++ -o test `pkg-config --cflags playerc++` test.cpp `pkg-config --libs playerc++`  -I/usr/local/include/player-2.1 -L/usr/lib */

#include <iostream>
#include <csignal>
#include <cstdlib> // for system()
#include <numeric>
#include <cstdio>

#include <iterator>
#include <fstream>

#include "libplayerc++/playerc++.h"

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "awi.h"
//#include "global.h"

#define OBSTACLE_HIT
// #define SPIKY_FWD_VEL
// #define SPIKY_TURN_VEL
// #define SMALL_L1
// #define SLOW_FWD_VEL
// #define AP_CONSTRAINTS
// #define BIG_DEL_VF


extern double gbest_x, gbest_y;
extern size_t ggen_idx, gind_idx, gbest_gen, gbest_ind, gbest_i;
//extern std::ofstream gclosest_ever_file;
extern FILE* gclosest_ever_fptr;

//extern std::ofstream sonar_file;

extern std::ofstream gsummary_file;
extern int gport;

std::ofstream gcons_file("cons.txt", std::ios::out);
std::ofstream gweights_file("w.txt", std::ios::out);
size_t gsim_idx;

double const tx = 12.0, ty = 12.0; // target position

void control(double * u, double * in, double * xreal);
void ann_control(double * u, double * in, double * xreal);

void show_pos2d(PlayerCc::PlayerClient & robot, 
		PlayerCc::Position2dProxy & pp) {
  robot.Read();
  std::cout << "pos2d :: " << pp.GetXPos() << "\t"
	    << pp.GetYPos() << "\t"
	    << pp.GetYaw() << "\t"
	    << pp.GetXSpeed() << "\t"
	    << pp.GetYSpeed() << "\t"
	    << pp.GetYawSpeed() << "\t"
	    << pp.GetStall() << std::endl;
}


void show_sonar(PlayerCc::PlayerClient & robot, 
		PlayerCc::SonarProxy & sp) {
  robot.Read();
  std::cout << "son :: ";
  for(size_t i = 0; i < 8; ++i) {
    std::cout << sp[i] << "\t";
  }
  std::cout << std::endl;
}


void show_pose(PlayerCc::PlayerClient & robot, 
	       PlayerCc::SimulationProxy & sim, char * name) {
  double xx = 0.0, yy = 0.0, aa = 0.0;
  robot.Read();
  sim.GetPose2d(name, xx, yy, aa);
  std::cout << "pose :: " << xx << "\t" << yy << "\t" << aa << std::endl;
}


void reset_stall(PlayerCc::PlayerClient & robot,
		 PlayerCc::Position2dProxy & pp) {
  while(pp.GetStall()) {
    //std::cout << "resetting stall ..." << std::endl;
    pp.SetSpeed(-0.1,0.0);
    robot.Read();
  }

  while((pp.GetXSpeed() != 0.0) || (pp.GetYSpeed() != 0.0)) {
    pp.SetSpeed(0.0, 0.0);
    robot.Read();
  }
}


void set_pose(double const x, double const y, double const a, 
	      PlayerCc::PlayerClient & robot, 
	      PlayerCc::Position2dProxy & pp, 
	      PlayerCc::SimulationProxy & sim, char * name) {
  //std::cout << "\nsetting pose ..." << std::endl;
  double xx = 0.0, yy = 0.0, aa = 0.0;
  robot.Read();
  sim.GetPose2d(name, xx, yy, aa);
  //reset_stall(robot, pp);
  while((xx != x) || (yy != y) || (aa != a)) {
    sim.SetPose2d(name, x, y, a);
    pp.SetOdometry(0.0, 0.0, a);
    pp.SetSpeed(0.0, 0.0);
    robot.Read();
    sim.GetPose2d(name, xx, yy, aa);
    //show_pose(robot, sim, name);
    //SHOW(pp.GetStall());
  }

  reset_stall(robot, pp);

  while(pp.GetXPos() != 0.0) {
    pp.SetOdometry(0.0, 0.0, a);
    robot.Read();
  }
  //std::cout << "done\n" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
void normalize(double * in, double * lower, 
	       double * upper, size_t const n){
  for(size_t i = 0; i < n; ++i){
    in[i] = 2*( (in[i] - lower[i]) / (upper[i] - lower[i]) 
		- 0.5);
  }
}

///////////////////////////////////////////////////////////////////////////////
struct abs_less_than{
  double val_;
  abs_less_than(double val):val_(val) {}
  bool operator()(double x) const { return fabs(x) < val_;}
};

///////////////////////////////////////////////////////////////////////////////
size_t num_sign_changes(std::vector<double> const & v){
  size_t ans = 0;
  for(size_t i = 1; i < v.size(); ++i){
    if(v[i]*v[i-1] < 0){
      ++ans;
    }
  }
  return ans;
}

///////////////////////////////////////////////////////////////////////////////
void fft(double * mag, double * data, size_t const n){
  gsl_fft_complex_radix2_forward(data, 1, n);
  for(size_t i = 0; i < n; ++i){
    mag[i] = hypot(data[2*i], data[2*i+1]);
  }

  double max_mag = *std::max_element(mag, mag+n);
  
  for(size_t i = 0; i < n; ++i){
    mag[i] /= max_mag;
  }
}


///////////////////////////////////////////////////////////////////////////////
bool all_almost_same(std::vector<double> & v, double tol){
  double mean = gsl_stats_mean(v.data(), 1, v.size());
  for(size_t i = 0; i < v.size(); ++i){
    v[i] -= mean;
  }
  
  return std::all_of(v.begin(), v.end(), abs_less_than(tol));
}

///////////////////////////////////////////////////////////////////////////////
size_t num_spikes_in(std::vector<double> const & v, double const mag){
  size_t count = 0;
  for(size_t i = 2; i < v.size(); ++i){
    if(fabs(v[i] - v[i-1]) > mag &&
       fabs(v[i-1] - v[i-2]) > mag && 
       (v[i] - v[i-1])*(v[i-1] - v[i-2]) < 0){
      ++count;
    }
  }
  return count;
}


///////////////////////////////////////////////////////////////////////////////
double calc_avoidance_penalty(std::vector<double> const & s, 
			      std::vector<double> const & prev_s) {
  std::vector<double>::const_iterator itr = 
    std::min_element(prev_s.begin()+1, prev_s.begin()+7);

  int i = static_cast<int>(std::distance(prev_s.begin()+1, itr));
				      
  // std::cout << i << std::endl;
  // for(size_t k = 1; k < 7; ++k) {
  //   std::cout << s[k] << "\t";
  // }
  // std::cout << std::endl;

  // for(size_t k = 1; k < 7; ++k) {
  //   std::cout << prev_s[k] << "\t";
  // }
  // std::cout << std::endl;

  if(s[i] > prev_s[i]) {
    //std::cout << "avoid +1" << std::endl;
    return 1.0/(0.1 + prev_s[i]);
  } else {
    //std::cout << "avoid -1" << std::endl;
    return -1.0/(0.1 + s[i]);
  }
}


///////////////////////////////////////////////////////////////////////////////
double calc_avoidance_penalty2(const double in[5], const double prev_in[5]) {
  size_t i = std::distance(prev_in, 
			   std::min_element(prev_in, prev_in+3));
  return (prev_in[i] - in[i])/(in[i] + 0.01);
}


///////////////////////////////////////////////////////////////////////////////
double calc_avoidance_penalty3(const double in[5], const double prev_in[5]) {
  double penalty = 0.0;
  for(size_t i = 0; i < 3; ++i) {
    penalty += (prev_in[i] - in[i])/(in[i] + 0.01);
  }
  return penalty;
}


///////////////////////////////////////////////////////////////////////////////
double calc_reaching_penalty(double const x, double const y, double const a,
			     double const px, double const py, double const pa) {
  double curr_dist = hypot(tx-x, ty-y);
  double prev_dist = hypot(tx-px, ty-py);
  double dist_penalty = 0.0;
  if(curr_dist < prev_dist) {
    dist_penalty = -0.5;
  } else {
    dist_penalty = 0.5;
  }

  double curr_ang_dist = (atan2(ty - y, tx - x) - a);
  double prev_ang_dist = (atan2(ty - py, tx - px) - pa);
  double ang_penalty = 0.0;
  if(curr_ang_dist < prev_ang_dist) {
    ang_penalty = -0.5;
  } else {
    ang_penalty = 0.5;
  }

  return (dist_penalty + ang_penalty)/(1+curr_dist);
}


double calc_reaching_penalty4(double * in, double * prev_in) {
  double p1 = 0.0;
  if(in[3] >= prev_in[3]){ // moved away
    p1 = 0.9;
  } else {
    p1 = -0.5;
  }

  double p2 = 0.0;
  if(in[4] >= prev_in[4]){ // moved away
    p2 = 0.9;
  } else {
    p2 = -0.5;
  }

  return p1+p2;
}


///////////////////////////////////////////////////////////////////////////////
double calc_avoidance_penalty4(const double in[5], const double prev_in[5],
			       double a, double prev_a) {
  double penalty = 0.0;
  for(size_t i = 0; i < 3; ++i) {
    penalty += (prev_in[i] - in[i])*exp(-in[i]);
  }

  // double p2 = 0.0;
  // size_t i = std::max_element(prev_in, prev_in+3);
  
  return penalty;
}


///////////////////////////////////////////////////////////////////////////////
double limit(double const x, double const min = -1.0, double const max = 1.0) {
  assert(min<max);
  return x<min?min:(x>max?max:x);
}


///////////////////////////////////////////////////////////////////////////////
void hack() {
  PlayerCc::PlayerClient robot("localhost", gport);
  PlayerCc::SimulationProxy sim(&robot);
  PlayerCc::Position2dProxy pp(&robot);
  PlayerCc::SonarProxy sp(&robot, 0);

  for(size_t i = 0; i < 10; i++) {
    robot.Read();
    pp.SetSpeed(-1.0, 0.0);
  }
}


///////////////////////////////////////////////////////////////////////////////
void cost_function(double * obj, double * constr, double *xreal, 
		   std::vector<std::vector<double> > & data, bool & hit, 
		   bool & reached, double & proximity, size_t & lifetime,
		   double & dist, double & lc, double & fc, double & rc) {
  
  //std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  std::cout << "\n==========> gen :: " << ggen_idx << "\t individual :: " 
	    << gind_idx << std::endl;

  // if(0 != gind_idx){
  //   return;
  // }

  PlayerCc::PlayerClient robot("localhost", gport);
  PlayerCc::SimulationProxy sim(&robot);
  PlayerCc::Position2dProxy pp(&robot);
  PlayerCc::SonarProxy sp(&robot, 0);

  double const x0 = -7.0;
  double const y0 = -7.0;
  double const a0 = PlayerCc::dtor(45);
  //double const a0 = PlayerCc::dtor(0);

  // control signal
  double u[2]={0.0, 0.0};
  double in[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double prev_in[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  prev_in[3] = hypot(tx - 0, ty - 0);
  prev_in[4] = (atan2((ty - y0), (tx - x0)) - a0);
  // double in_upper[5] = {5.0, 5.0, 5.0, 17.0, M_PI};
  // double in_lower[5] = {0.0, 0.0, 0.0, 0.0, -M_PI};
  double in_upper[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  double in_lower[5] = {0.0, 0.0, 0.0, 0.0, -1.0};

  double x = 0.0, y = 0.0, a = 0.0; 
  double prev_x = 0.0, prev_y = 0.0, prev_a = 0.0; 
  std::vector<double> s(16, 0.0);

  double obstacle = 0.0; // obstacle avoidance
  double target = 0.0;   // target nearness 
  
  reached = false;
  hit = false;
  dist = 0.0;

  char robname[] = "robot1";

  set_pose(x0, y0, a0, robot, pp, sim, robname);
  //show_pos2d(robot, pp);
  //show_sonar(robot, sp);

  size_t const tmax = 500; // robot lifetime
  size_t tfinal = tmax-1;

  double reaching_penalty = 0.0;
  double avoidance_penalty = 0.0;

  double target_penalty = 0.0;
  double obstacle_penalty = 0.0;

  // BEGIN SIMULATION 
  for(size_t i = 0; i < tmax; ++i) {
    gsim_idx = i;
    robot.Read();

    if(pp.GetStall()) {
      std::cout << "obstacle hit at " << i << "..." << std::endl;
      hit = true;
      tfinal = i;
      //std::cout << x << "\t" << y << std::endl;

      // auto itr = std::min_element(s.begin(), s.end());
      // std::cout << "min sonar :: " << std::distance(s.begin(), itr) 
      // 		<< "\t" << "value :: " << *itr << std::endl;
      // show(s);

      //show_sonar(robot, sp);

      break;
    }

    a = pp.GetYaw(); x = pp.GetXPos(); y = pp.GetYPos();
    data.push_back(std::vector<double>());
    data[i].push_back(x);
    data[i].push_back(y);
    data[i].push_back(a);

    for(size_t k = 0; k < 16; ++k) { s[k] = sp[k]; }

    // {
    //   for(size_t kk = 0; kk < 8; ++kk) {
    // 	sonar_file << s[kk] << "\t";
    //   }
    //   sonar_file << "\n";
    // }

    dist += hypot(x - prev_x, y - prev_y);
    //std::cout << "i :: " << i << "\t" << "dist :: " << dist << std::endl;

    prev_x = x; prev_y = y; 

    // //in[0] = std::min(s[1], s[2]); 
    // in[0] = *std::min_element(s.begin(), s.begin()+3); // left
    // in[1] = std::min(s[3], s[4]); 
    // //in[2] = std::min(s[5], s[6]); 
    // in[2] = *std::min_element(s.begin()+5, s.begin()+8); // right
    // in[3] = hypot(tx - x, ty - y); // distance to target
    // in[4] = (atan2((ty - y), (tx - x)) - a);
    // if(in[4] > +M_PI) {in[4] -= 2.0*M_PI;}
    // if(in[4] < -M_PI) {in[4] += 2.0*M_PI;}

    in[0] = ((s[0] + s[1] + s[2]) / 15.0) ; // left clearance
    in[1] = ((s[3] + s[4]) / 10.0); 
    in[2] = ((s[5] + s[6] + s[7]) / 15.0) ; // right clearance
    in[3] = (hypot(tx - x, ty - y) / 17.0); // distance to target
    in[4] = (atan2((ty - y), (tx - x)) - a);
    if(in[4] > +M_PI) {in[4] -= 2.0*M_PI;}
    if(in[4] < -M_PI) {in[4] += 2.0*M_PI;}
    in[4] /= M_PI;

    {
      for(size_t kk = 0; kk < 5; ++kk) {
	data[i].push_back(in[kk]);
      }

      lc += (in[0]);
      fc += (in[1]);
      rc += (in[2]);
    }

    obstacle  += ( ( exp(-s[0]) + exp(-s[1]) + exp(-s[2]) + exp(-s[3]) + 
		     exp(-s[4]) + exp(-s[5]) + exp(-s[6]) + exp(-s[7]) ) / 8.0 );

    // target nearness
    target += ((in[3] + fabs(in[4])) / 2.0);

    //if((hypot(tx - x, ty - y) < 0.1) && (fabs(in[4]) < 0.25)) {
    if((hypot(tx - x, ty - y) < 0.1)) {
      std::cout << "i :: " << i << std::endl;
      std::cout << "\033[0;31m reached ... \033[0m" << std::endl;
      std::cout << "ggen_idx " << ggen_idx << std::endl;
      std::cout << "gind_idx " << gind_idx << std::endl;
      reached = true;
      tfinal = i;
      break;
    }

//     if(pp.GetStall()) {
//       std::cout << "obstacle hit at " << i << "..." << std::endl;
//       hit = true;
//       tfinal = i;
//       std::cout << x << "\t" << y << std::endl;
//       // auto itr = std::min_element(s.begin(), s.end());
//       // std::cout << "min sonar :: " << std::distance(s.begin(), itr) 
//       // 		<< "\t" << "value :: " << *itr << std::endl;
//       // show(s);

//       break;
//     }

    normalize(in, in_lower, in_upper, 5);

    //ann_control(u, in, xreal);
    control(u, in, xreal);
    
    u[0] = limit(u[0], -1.0, 1.0); // forward
    u[1] = limit(u[1], -1.0, 1.0); // rotational
    
    //u[0] = -1; u[1] = 0.0;
    pp.SetSpeed(u[0], u[1]);
    //pp.SetSpeed(0.1, 0.0);
    //std::cout << u[0] << "\t" << u[1] << "\n";

    data[i].push_back(u[0]);
    data[i].push_back(u[1]);

    {
      for(size_t kk = 0; kk < 8; ++kk) {
	data[i].push_back(s[kk]);
      }
    }

  } // END OF SIMULATION 
  
  lifetime = tfinal;
  proximity = hypot(tx - x, ty - y);

  double hit_penalty = 1.0;
  if(hit) {
    hit_penalty = 2.0;
  }

  if(!reached) {
    constr[0] = -1.0 * hit_penalty * ( (hypot(tx - x, ty - y)/ 17.0) +
  	(fabs((atan2((ty - y), (tx - x)) - a)) / M_PI) );
  } else { 
    constr[0] = 0.0;
  }

  // maximize sonar readings i.e. stay away from obstacles
  //obj[0] = -(obstacle_farness); 
  //obj[0] = -(obstacle_farness - del_v_sum); 
  //obj[0] = avoidance_penalty;
  obj[0] = (obstacle)/double(tfinal); 
  //obj[0] = (obstacle_penalty) / double(tfinal); 

  // minimize target farness
  //obj[1] = (target_farness);
  obj[1] = (dist)/double(tfinal);
  //obj[1] = (target_farness + del_v_sum);
  //obj[1] = reaching_penalty;
  //obj[1] = (awayness)/double(tfinal);
  //obj[1] = (target_penalty) / double(tfinal);

  // SHOW(obj[0]);
  // SHOW(obj[1]);
  // SHOW(constr[0]);

  //sonar_file << "\n";

  return ;
}
