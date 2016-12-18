#include <gsl/gsl_rng.h>
#include <iostream>
#include <fstream>

double gbest_x, gbest_y;
int ggen_idx, gind_idx, gbest_gen, gbest_ind;
FILE* gclosest_ever_fptr;
std::ofstream gsummary_file("summary.txt", std::ios::out);
int gport = 6665;

void cost_function(double * obj, double * constr, double *xreal);

int fuzzy_tester(){
  double obj[2];
  double constr[1];
  double xreal[2916];

  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc(T);

  gsl_rng_env_setup();
  
  for(size_t i = 0; i < 2916; ++i){
    xreal[i] = -1.0 + 2.0*gsl_rng_uniform(r);
  }

  // for(size_t i = 0; i < 2916; ++i){
  //   std::cout << xreal[i] << std::endl;
  // }

  cost_function(obj, constr, xreal);
  return 0;
}

int ann_tester(){
  size_t const nvars = 34;
  double obj[2];
  double constr[1];
  double xreal[nvars];

  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc(T);

  gsl_rng_env_setup();
  
  for(size_t i = 0; i < nvars; ++i){
    xreal[i] = -1.0 + 2.0*gsl_rng_uniform(r);
  }

  cost_function(obj, constr, xreal);
  return 0;
}


int main(int argc, char ** argv){
  gport = atoi(argv[2]);
  ann_tester();
}
  
