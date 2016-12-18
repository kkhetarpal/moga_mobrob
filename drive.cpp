#include "libplayerc++/playerc++.h"

#include <random>

int main() {

  std::mt19937 gen();
  //std::uniform_real_distribution<> dis(0,1);

  std::stringstream ss;
  std::string command =  "player -p 6665 simple.cfg > /dev/null 2>&1 &";
  system(command.c_str());
  system("sleep 1");
  //system("sleep 10"); // extra time to enable trails
  
  PlayerCc::PlayerClient robot("localhost", 6665);
  PlayerCc::SimulationProxy sim(&robot);
  PlayerCc::Position2dProxy pp(&robot);
  PlayerCc::SonarProxy sp(&robot, 0);

  double uf = 0.1, ut = 0.01;
  double a = 0.0;
  double lc = 0.0, fc = 0.0, rc = 0.0;

  std::cout.setf(std::ios::showpos | std::ios::scientific | std::ios::showpoint);
  std::cout << uf << ut << std::endl;

  for(size_t i = 0; i < 500; i++) {
    robot.Read();

    // lc = (sp[0] + sp[1] + sp[2])/15.0;
    // rc = (sp[5] + sp[6] + sp[7])/15.0;

    // std::cout << lc << "\t" << rc << std::endl;

    if(pp.GetStall()) {
      std::cout << "obstacle hit ..." << std::endl;
      break;
    }

    pp.SetSpeed(0.5, -0.1);
  }

  {
    std::cout << pp.GetStall() << std::endl;
    
  }
}
