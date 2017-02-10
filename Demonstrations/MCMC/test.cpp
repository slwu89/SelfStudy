#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(cpp11)]]
#include <progress.hpp>
#include <random>

using namespace Rcpp;


// MCMC utility functions
// [[Rcpp::export]]
arma::vec mvrnorm_samp(arma::vec mu, arma::mat sigma) {
  arma::rowvec Y = rnorm(sigma.n_cols,0,1);
  arma::rowvec out = mu.t() + Y * arma::chol(sigma);
  return(out.t());
}


// std version of RNG
std::mt19937 engine(42); //set seed of uniform RNG
std::normal_distribution<> normal(0.0,1.0);

// [[Rcpp::export]]
arma::vec boostTest(arma::vec mu, arma::mat sigma) {
  arma::rowvec Y(sigma.n_cols);
  for(int i = 0; i < sigma.n_cols; i++){
    Y(i) = normal(engine);
  }
  arma::rowvec out = mu.t() + Y * arma::chol(sigma);
  return(out.t());
}