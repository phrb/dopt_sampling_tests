#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_biased_sample(NumericMatrix cdf_data, int samples, int factors){
  NumericMatrix output(samples, factors);
  float runif_sample;
  int k;

  for(int i = 0; i < output.nrow(); i++){
    for(int j = 0; j < output.ncol(); j++){
      k = 0;
      runif_sample = runif(1)[0];

      while(k < cdf_data.nrow()){
        if(cdf_data(k, 1) >= runif_sample){
          output(i, j) = cdf_data(k, 0);
          break;
        }

        k++;
      }
    }
  }
  return output;
}
