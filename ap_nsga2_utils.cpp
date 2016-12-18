# include "global.h"
#include <iostream>

size_t num_reaches(population * pop){
  size_t ans = 0;
  for(size_t i = 0; i < popsize; ++i){
    //std::cout << "constr_violation :: " << pop->ind[i].constr_violation << std::endl;
    if(pop->ind[i].reached){
      ans++;
    }
  }
  return ans;
}


size_t num_hits(population * pop){
  size_t ans = 0;
  for(size_t i = 0; i < popsize; ++i){
    //std::cout << "constr_violation :: " << pop->ind[i].constr_violation << std::endl;
    if(pop->ind[i].hit){
      ans++;
    }
  }
  return ans;
}
    
