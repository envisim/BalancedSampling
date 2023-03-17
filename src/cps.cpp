#include <Rcpp.h>
#include "cps-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.cps_cpp)]]
Rcpp::IntegerVector cps_cpp(
  int cpsMethod,
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  int p = x.nrow();

  if (prob.length() != N)
    std::invalid_argument("prob an x does not match");

  Cps cps(REAL(prob), REAL(x), nullptr, N, p, intToCpsMethod(cpsMethod), bucketSize, method, eps);
  cps.run();

  Rcpp::IntegerVector sample(cps.sample, cps.sample + cps.sampleSize);

  return sample;
}

// [[Rcpp::export(.cps_random_cpp)]]
Rcpp::IntegerVector cps_random_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::NumericVector &random,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  int p = x.nrow();

  if (prob.length() != N)
    std::invalid_argument("prob an x does not match");
  if (random.length() != N)
    std::invalid_argument("random an x does not match");

  Cps cps(
    REAL(prob),
    REAL(x),
    REAL(random),
    N,
    p,
    CpsMethod::scpscoord,
    bucketSize,
    method,
    eps
  );
  cps.run();

  Rcpp::IntegerVector sample(cps.sample, cps.sample + cps.sampleSize);

  return sample;
}
