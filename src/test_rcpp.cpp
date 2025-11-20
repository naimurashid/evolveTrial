// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Test Rcpp setup
//'
//' Simple test function to verify Rcpp and RcppArmadillo are working.
//'
//' @param x A numeric vector
//' @return Sum of the vector
//' @export
// [[Rcpp::export]]
double test_rcpp_sum(const arma::vec& x) {
  return arma::sum(x);
}

//' Test Armadillo matrix operations
//'
//' @param n Matrix dimension
//' @return Mean of random normal matrix
//' @export
// [[Rcpp::export]]
double test_armadillo_random(int n) {
  arma::mat A = arma::randn<arma::mat>(n, n);
  return arma::mean(arma::mean(A));
}
