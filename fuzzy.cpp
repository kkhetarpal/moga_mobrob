#include <cstddef>
#include <vector>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

//#include "ap_utils.h"
#include "awi.h"

extern std::ofstream gcons_file;
extern std::ofstream gweights_file;
extern size_t gsim_idx;


typedef struct {
  size_t ni; // number of inputs
  size_t no; // number of outputs
  size_t nr; // number of rules
} info_t;

double ap_gsl_prod(const gsl_vector *v){
  double prod = 1.0;
  for(size_t i = 0; i < v->size; ++i) prod *= gsl_vector_get(v, i);
  return prod;
}

double ap_gsl_sum(const gsl_vector *v){
  double sum = 0.0;
  for(size_t i = 0; i < v->size; ++i) sum += gsl_vector_get(v, i);
  return sum;
}


///////////////////////////////////////////////////////////////////////////////
void memberships(gsl_matrix * mu, gsl_matrix const * x, gsl_matrix const * c,
                 gsl_matrix const * s){
  size_t const np = x->size1;
  size_t const ni = x->size2;
  size_t const nr = c->size1;
  
  for(size_t i = 0; i < np; ++i){
    for(size_t j = 0; j < nr; ++j){
      for(size_t k = 0; k < ni; ++k){
        double const xik = gsl_matrix_get(x, i, k);
        double const cjk = gsl_matrix_get(c, j, k);
        double const sjk = gsl_matrix_get(s, j, k);
        double z = (xik-cjk)/sjk; 
        gsl_matrix_set(mu, i*nr+j, k, exp(-0.5*z*z));
      }
    }
  }
  return;
}


///////////////////////////////////////////////////////////////////////////////
void truth(gsl_matrix *tau, const gsl_matrix *mu, size_t const np,
           size_t const nr){
  for(size_t p = 0; p < np; ++p){
    for(size_t r = 0; r < nr; ++r){
      gsl_vector_const_view mu_row_view =
        gsl_matrix_const_row(mu, p*nr+r);
      double value = ap_gsl_prod(&mu_row_view.vector);
      gsl_matrix_set(tau, p, r, value);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
void weights(gsl_matrix *w, const gsl_matrix *tau){
  // divide each row by the sum of row elements
  gsl_matrix_memcpy(w, tau); // w <-- tau
  for(size_t r = 0; r < w->size1; ++r){
    gsl_vector_view v = gsl_matrix_row(w, r);
    double scale = 1.0/ap_gsl_sum(&v.vector);
    gsl_vector_scale(&v.vector, scale);
  }

  // {
  //   for(size_t r = 0; r < w->size1; ++r){
  //     double sum = 0.0;
  //     for(size_t c = 0; c < w->size2; ++c){
  // 	sum += gsl_matrix_get(w, r, c);
  //     }
  //     std::cout << "sum :: " << sum << std::endl;
  //   }
  // }

  // {
  //   gweights_file << gsim_idx << "\t";
  //   for(size_t ii = 0; ii < w->size1; ++ii){
  //     for(size_t jj = 0; jj < w->size2; ++jj){
  // 	gweights_file << gsl_matrix_get(w, ii, jj) << "\t";
  //     }
  //     gweights_file << std::endl;
  //   }
  // }
}


///////////////////////////////////////////////////////////////////////////////
void fuzzy_outputs(gsl_matrix *y, const gsl_matrix *x, const gsl_matrix *w,
                   const gsl_matrix *cons, info_t const info){
  // {
  //   std::cout << "x :: " << "\t";
  //   double one = 1.0;
  //   std::cout << one << "\t";
  //   for(size_t ii = 0; ii < x->size2; ++ii){
  //     std::cout << gsl_matrix_get(x, 0, ii) << "\t";
  //   }
  //   std::cout << std::endl;
  // }

  gsl_matrix *H = gsl_matrix_alloc(info.no, info.ni+1);
  gsl_matrix *Hr = gsl_matrix_alloc(info.no, info.ni+1);
  gsl_vector *z = gsl_vector_alloc(1+info.ni);
  gsl_vector_set(z, 0, 1.0);    /* bias term */
  for(size_t i = 0; i < x->size1; ++i){
    gsl_matrix_set_all(H, 0.0);
    for(size_t r = 0; r < info.nr; ++r){
      gsl_matrix_const_view cons_view =
        gsl_matrix_const_submatrix(cons, r*info.no, 0, info.no, info.ni+1);
      gsl_matrix_memcpy(Hr, &cons_view.matrix);
      gsl_matrix_scale(Hr, gsl_matrix_get(w,i,r));
      gsl_matrix_add(H, Hr);  /* H += Hr */
    }
    
    // {
    //   std::cout << "H :: " << "\t";
    //   for(size_t ii = 0; ii < 1; ++ii){
    // 	for(size_t jj = 0; jj < H->size2; ++jj){
    // 	  std::cout << gsl_matrix_get(H, ii, jj) << "\t";
    // 	  gcons_file << gsl_matrix_get(H, ii, jj) << "\t";
    // 	}
    // 	std::cout << std::endl;
    // 	gcons_file << std::endl;
    //   }
    // }
    
    {
      gsl_vector_view zsub_view = gsl_vector_subvector(z, 1, info.ni);
      gsl_vector_const_view xrow_view = gsl_matrix_const_row(x, i);
      gsl_vector_memcpy(&zsub_view.vector, &xrow_view.vector);
    }
    gsl_vector_view yrow = gsl_matrix_row(y, i);
    gsl_blas_dgemv(CblasNoTrans, 1.0, H, z, 0.0, &(yrow.vector));
  }

  //show(H->data, 2, 6);

  gsl_matrix_free(H);
  gsl_matrix_free(Hr);
  gsl_vector_free(z);
}


///////////////////////////////////////////////////////////////////////////////
void tsk(gsl_matrix * y, const gsl_matrix * x, const gsl_matrix * centers,
         const gsl_matrix * sigmas, const gsl_matrix * consequents, info_t
         const info){
  size_t const np = x->size1;
  size_t const ni = x->size2;
  size_t const nr = centers->size1;
  auto mu = gsl_matrix_alloc(np*nr, ni);
  memberships(mu, x, centers, sigmas);
  gsl_matrix *tau = gsl_matrix_alloc(x->size1, nr);
  truth(tau, mu, np, nr);
  gsl_matrix * w = gsl_matrix_alloc(x->size1, nr);
  weights(w, tau);
  fuzzy_outputs(y, x, w, consequents, info);
  gsl_matrix_free(mu);
  gsl_matrix_free(tau);
  gsl_matrix_free(w);
}


///////////////////////////////////////////////////////////////////////////////
void combinations(std::vector<std::vector<double> > const & sets,
                  gsl_matrix * combs){
  size_t const nrows = sets.size();

  size_t ncombs = 1;
  //for(auto i: sets){ncombs *= i.size();}
  for(size_t i = 0; i < sets.size(); ++i){
    ncombs *= sets[i].size();
  }
                    
  std::vector<size_t> count(nrows, 0);
  size_t i = 0;
  for(; i < ncombs;){
    for(size_t j = 0; j < nrows; ++j){
      gsl_matrix_set(combs, i, j, sets[j][count[j]]);
    }
    ++i;

    count[nrows-1] = i % sets[nrows-1].size(); // increment right most
    
    if(0 == count[nrows-1]){ // if right most index wraps back to 0
      for(size_t k = 1; k < nrows; ++k){ // traverse count from right to left
        if(0 == count[nrows-k]){
          count[nrows-1-k] = (count[nrows-1-k] + 1) % sets[nrows-1-k].size();
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void control(double * u, double * in, double * xreal){
  // grng = gsl_rng_alloc(gsl_rng_default);
  // gsl_rng_set(grng, time(NULL));

  size_t const ni = 5;
  size_t const no = 2;
  info_t info{ni, no, 0};

  std::vector<size_t> nzones{3, 3, 3, 3, 3};
  std::vector<std::vector<double> > c(ni);
  std::vector<std::vector<double> > s(ni);

  for(std::size_t i = 0; i < ni; ++i){ 
    for(std::size_t j = 0; j < nzones[i]; ++j){
      c[i].push_back(-1.0 + j*2.0/(nzones[i]-1));
    }
  }

  {
    double const overlap = 0.2;
    double const konst = sqrt(-8*log(overlap));
    for(std::size_t i = 0; i < ni; ++i){
      double const sig = (c[i][1] - c[i][0]) / konst;
      for(std::size_t j = 0; j < nzones[i]; ++j){
        s[i].push_back(sig);
      }
    }
  }

  //display(c);
  //display(s);

  // size_t const nr = std::accumulate(begin(nzones), end(nzones),
  //                                   1, std::multiplies<size_t>());
  
  size_t nr = 1;
  for(size_t i = 0; i < ni; ++i){
    nr *= nzones[i];
  }

  info.nr = nr;
  //std::cout << "nr :: " << nr << std::endl;
  //SHOW(nr*no*(ni+1));

  gsl_matrix * centers = gsl_matrix_alloc(nr, ni);
  gsl_matrix * sigmas = gsl_matrix_alloc(nr, ni);

  combinations(c, centers);
  combinations(s, sigmas);

  gsl_matrix * cons = gsl_matrix_alloc(no*nr, ni+1);
  //gsl_matrix_set_all(cons, 1.0);
  for(size_t i = 0, k = 0; i < no*nr; ++i){
    for(size_t j = 0; j < (ni+1); ++j, ++k){
      gsl_matrix_set(cons, i, j, xreal[k]);
      //std::cout << i << "\t" << j << "\t" << k << std::endl;
    }
  }

  // for(size_t i = 0, k = 0; i < no*nr; ++i){
  //   for(size_t j = 0; j < (ni+1); ++j, ++k){
  //     std::cout << gsl_matrix_get(cons, i, j) << "\t";
  //   }
  //   std::cout << std::endl;
  // }

  gsl_matrix * input = gsl_matrix_alloc(1, ni);
  gsl_matrix * output = gsl_matrix_alloc(1, no);

  //gsl_matrix_set_all(input, 1.0);
  for(size_t i = 0; i < ni; ++i){
    gsl_matrix_set(input, 0, i, in[i]);
  } 

  tsk(output, input, centers, sigmas, cons, info);
  
  for(size_t i = 0; i < no; ++i){
    //std::cout << gsl_matrix_get(output, 0, i) << std::endl;
    u[i] = gsl_matrix_get(output, 0, i);
  }

  // std::cout << "u :: ";
  // for(size_t i = 0; i < no; ++i){
  //   std::cout << u[i] << "\t";
  // }
  // std::cout << std::endl;

  //AP_GSL_SHOW(output);

  //gsl_rng_free(grng);
  gsl_matrix_free(centers);
  gsl_matrix_free(sigmas);
  gsl_matrix_free(cons);
  gsl_matrix_free(input);
  gsl_matrix_free(output);

  // u[0] = 0.2;
  // u[1] = 0.0;
  return;
}
