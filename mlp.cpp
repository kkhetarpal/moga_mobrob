#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

//#include "ap_utils.h"

//extern std::vector<double> gxtrn;
//extern std::vector<double> gytrn;
//extern size_t gxdim;
//extern size_t gydim;

///////////////////////////////////////////////////////////////////////////////
std::vector<double> ann(std::vector<double> const & x,
                        std::vector<size_t> const & nn,
                        std::vector<std::function<double(double)> >const& act,
                        std::vector<double> const & vars){      
  std::vector<double> z(x);
  size_t const nlayers = nn.size();

  std::vector<double> y;
  for(size_t i = 0; i < (nlayers-1); ++i){
    //SHOW(i);
    z.push_back(1.0); // add bias term
    gsl_vector_view in = gsl_vector_view_array(z.data(), z.size());
    //AP_GSL_SHOW(&in.vector);
    
    size_t offset = 0;
    for(size_t j = 0; j < i; ++j){offset += nn[j+1]*(nn[j]+1);} 
    gsl_matrix_const_view w =
      gsl_matrix_const_view_array(vars.data()+offset, nn[i+1], nn[i]+1);
    
    //ap_gsl_show(&w.matrix);
    //z.resize(nn[i+1]);
    y.resize(nn[i+1], 0.0);
    gsl_vector_view out = gsl_vector_view_array(y.data(), y.size());
    //AP_GSL_SHOW(&out.vector);
    gsl_blas_dgemv(CblasNoTrans, 1, &w.matrix, &in.vector, 0, &out.vector);
    //AP_GSL_SHOW(&out.vector);
    //display(y);
    std::transform(y.begin(), y.end(), y.begin(), act[i]);
    //display(y);
    z = y;
  }
  return y;
}


// ///////////////////////////////////////////////////////////////////////////////
// double mismatch(double * xreal, std::vector<double> const & xtrn,
//                 std::vector<double> const & ytrn){
//   auto sigmoid = [](double xx){return 1/(1+exp(-xx));};
//   auto unity = [](double xx){return xx;};
//   std::vector<std::function<double(double)> > act{sigmoid, unity};
//   std::vector<double> y(1);
//   std::vector<double> yd(1);
//   std::vector<double> x(1);
//   std::vector<double> vars(xreal, xreal+16);

//   double e = 0;
//   for(size_t i = 0; i < 50; ++i){
//     std::copy(std::begin(xtrn) + i*1, std::begin(xtrn) + (i+1)*1,
//               std::begin(x));
//     std::copy(std::begin(ytrn) + i*1, std::begin(ytrn) + (i+1)*1,
//               std::begin(yd));
//     y = ann(x, {1,5,1}, act, vars);
//     e += one_norm(subtract(y, yd));
//   }
//   return e/50.0;
// }

// void doit(double * obj, double * xreal){
//   obj[0] = mismatch(xreal, gxtrn, gytrn);
  
//   obj[1] = 0;
//   for(size_t i = 0; i < 16; ++i){
//     obj[1] += fabs(xreal[i]);
//   }
// }


double sigmoid(double x){
  return 1.0/(1.0+exp(-x));
}

double unity(double x){
  return x;
}

double mysig(double const x) {
  return 1.7159*tanh(2.0*x/3.0);
}


size_t calc_nvars(std::vector<size_t> const & nn) {
  size_t nvars = 0;
  for(size_t i = 0; i < (nn.size()-1); ++i) {
    nvars += nn[i+1]*(nn[i]+1);
  }
  return nvars;
}


void ann_control(double * u, double * in, double * xreal){
  std::vector<double> x(in, in+5);

  // std::vector<size_t> nn{5,4,2}; 5 inputs - 4 neurons in hidden
  // layer - 2 outputs

  // std::vector<size_t> nn{5,4,2};
  // std::vector<double> vars(xreal, xreal+34);

  std::vector<size_t> nn{5,10,2};
  size_t const nvars = calc_nvars(nn);
  std::vector<double> vars(xreal, xreal+nvars);

  // activation function 
  std::vector<std::function<double(double)> > 
    act{mysig, mysig, mysig, unity};
  //act{sigmoid, sigmoid, unity};

  std::vector<double> uu = ann(x, nn, act, vars);
  u[0] = uu[0];
  u[1] = uu[1];
}
