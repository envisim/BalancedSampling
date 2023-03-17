#include <stdexcept>
#include <Rcpp.h>

#include "CubeStratifiedClass.h"

// [[Rcpp::export(.cube_stratified_cpp)]]
Rcpp::IntegerVector cube_stratified_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &strata,
  const double eps
) {
  int N = x.nrow();
  int p = x.ncol();

  if (prob.length() != N)
    std::invalid_argument("prob and x does not match");
  if (N != strata.length())
    std::range_error("strata and x does not match");

  CubeStratified cube(
    REAL(prob),
    REAL(x),
    INTEGER(strata),
    N,
    p,
    eps
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}

// [[Rcpp::export(.lcube_stratified_cpp)]]
Rcpp::IntegerVector lcube_stratified_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &xbalance,
  Rcpp::NumericMatrix &xspread,
  Rcpp::IntegerVector &strata,
  const int bucketSize,
  const int method,
  const double eps
) {
  int N = xbalance.nrow();
  int p = xbalance.ncol();
  int pxs = xspread.nrow();

  if (prob.length() != N)
    std::invalid_argument("prob and x does not match");
  if (N != strata.length())
    std::range_error("strata and x does not match");
  if (N != xspread.length())
    std::range_error("xspread and xbal does not match");

  CubeStratified cube(
    REAL(prob),
    REAL(xbalance),
    REAL(xspread),
    INTEGER(strata),
    N,
    p,
    pxs,
    eps,
    bucketSize,
    method
  );

  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}
