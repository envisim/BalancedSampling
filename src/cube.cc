#include <Rcpp.h>
#include "CubeClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.cube_cpp)]]
Rcpp::IntegerVector cube_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  double eps
) {
  int N = x.nrow();
  int p = x.ncol();

  if (prob.length() != N)
    std::invalid_argument("prob and x does not match");

  Cube cube(REAL(prob), REAL(x), N, p, eps);
  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}

// [[Rcpp::export(.lcube_cpp)]]
Rcpp::IntegerVector lcube_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &xbal,
  Rcpp::NumericMatrix &xspread,
  int bucketSize,
  int method,
  double eps
) {
  int N = xbal.nrow();
  int pbal = xbal.ncol();
  int pspread = xspread.nrow();

  if (N != xspread.ncol())
    std::invalid_argument("xbal and xspread does not match");
  if (prob.length() != N)
    std::invalid_argument("prob and x does not match");

  Cube cube(
    REAL(prob),
    REAL(xbal),
    N,
    pbal,
    eps,
    REAL(xspread),
    pspread,
    bucketSize,
    method
  );
  cube.Run();

  Rcpp::IntegerVector sample(cube.sample.begin(), cube.sample.end());

  return sample;
}
