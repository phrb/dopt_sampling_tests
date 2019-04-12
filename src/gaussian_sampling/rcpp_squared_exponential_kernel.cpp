#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_squared_exponential_kernel(NumericMatrix x,
                                              NumericMatrix y,
                                              float amplitude,
                                              NumericVector lengthscale){
  NumericMatrix output(x.nrow(), y.nrow());
  float distance;

  for(int i = 0; i < x.nrow(); i++) {
    for(int j = 0; j < y.nrow(); j++) {
      distance = 0;

      for(int k = 0; k < x.ncol(); k++) {
        distance += (((x(i, k) - y(j, k)) * (x(i, k) - y(j, k))) /
                      (lengthscale(k) * lengthscale(k)));
      }

      output(i, j) = amplitude * amplitude * exp(-0.5 * distance);
    }
  }
  return(output);
}
