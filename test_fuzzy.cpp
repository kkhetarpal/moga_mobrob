#include <cstddef>
#include <gsl/gsl_rng.h>

#include "awi.h"

const gsl_rng_type* grngType;            
gsl_rng* grng;                             

void control(double * u, double * in, double * xreal);

void all_weights_zero() {
  double u[2] = {-1.0, -1.0};
  double in[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double xreal[2916];

  // if no = 2, ni = 5, nr = 243 then nw = nr*no*(ni+1) = 2916
  for(size_t i = 0; i < 2916; ++i) { 
    xreal[i] = 0.0; 
  };

  for(size_t i = 0; i < 10; ++i) { 
    for(size_t j = 0; j < 5; ++j) {
      in[j] = -1.0 + 2.0*gsl_rng_uniform(grng); 
    }
    control(u, in, xreal);
    show(u, 2);
  }
}

void final_rule_same_if_all_rules_same() {

  double u[2] = {-1.0, -1.0};
  double in[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double xreal[2916];
  
  double rule[2*(5+1)];
  for(size_t i = 0; i < 12; ++i) {
    rule[i] = -1.0 + 2.0*gsl_rng_uniform(grng); 
  }

  size_t const nr = 243;
  for(size_t i = 1; i < nr; ++i) {
    for(size_t j = 0; j < 12; ++j) {
      xreal[i*12+j] = rule[j];
    }
  }
  
  show(rule, 2, 6);

  // turn on display in fuzzy_output()

  control(u, in, xreal);

}


int main() {
  gsl_rng_env_setup();
  grngType = gsl_rng_default;
  grng = gsl_rng_alloc(grngType);
  const unsigned long seed = time(NULL);
  gsl_rng_set(grng,  seed);

  //all_weights_zero();
  
  //final_rule_same_if_all_rules_same();

  gsl_rng_free(grng);
  return 0;
}
  
