// mcmc.cpp
#include <Rcpp.h>
using namespace Rcpp;

// variables
int it, nSim;
long double param_alpha;


// [[Rcpp::export(name = "mcmc")]]
void test(long double current_alpha, int nSim){
  Rprintf("%s %s %s %s %s \n", "nSim", "alpha", "gamma", "Gta", "Gtb");
  for(int it=0; it<nSim; it++){
    param_alpha = current_alpha + R::rnorm(0.0, 0.001);
    while(param_alpha < 0){
      param_alpha = current_alpha + R::rnorm(0.0, 0.001);
    }
  }
  Rprintf("%d %.10Lf\n", nSim, param_alpha);
}