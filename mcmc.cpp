// mcmc.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = "mcmc")]]
double test(NumericVector v){
  Rprintf("%s %s %s %s %s \n", "nSim", "alpha", "gamma", "Gta", "Gtb");
  double sum = 0;
  for(int i=0; i<v.length(); ++i){
    sum += v[i];
  }
  return(sum);
}