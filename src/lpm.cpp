#include <Rcpp.h>
#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.lpm_cpp)]]
Rcpp::IntegerVector lpm_cpp(
  int lpMethod,
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

  Lpm lpm(REAL(prob), REAL(x), N, p, intToLpmMethod(lpMethod), bucketSize, method, eps);
  lpm.run();

  Rcpp::IntegerVector sample(lpm.sample, lpm.sample + lpm.sampleSize);

  return sample;
}

// [[Rcpp::export(.lpm_int_cpp)]]
Rcpp::IntegerVector lpm_int_cpp(
  int lpMethod,
  int n,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method
) {
  int N = x.ncol();
  int p = x.nrow();

  Lpm lpm(n, REAL(x), N, p, intToLpmMethod(lpMethod), bucketSize, method);
  lpm.run();

  Rcpp::IntegerVector sample(lpm.sample, lpm.sample + lpm.sampleSize);

  return sample;
}

// [[Rcpp::export(.rpm_cpp)]]
Rcpp::IntegerVector rpm_cpp(
  Rcpp::NumericVector &prob,
  double eps
) {
  int N = prob.length();

  IndexList* idx = new IndexList(N);
  double* probabilities = new double[N];

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  Lpm lpm(probabilities, nullptr, idx, N, LpmMethod::rpm, eps);
  lpm.run();

  Rcpp::IntegerVector sample(lpm.sample, lpm.sample + lpm.sampleSize);

  delete idx;
  delete[] probabilities;

  return sample;
}

// [[Rcpp::export(.spm_cpp)]]
Rcpp::IntegerVector spm_cpp(
  Rcpp::NumericVector &prob,
  double eps
) {
  int N = prob.length();

  IndexList* idx = new IndexList(N);
  double* probabilities = new double[N];

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  Lpm lpm(probabilities, nullptr, idx, N, LpmMethod::spm, eps);
  lpm.run();

  Rcpp::IntegerVector sample(lpm.sample, lpm.sample + lpm.sampleSize);

  delete idx;
  delete[] probabilities;

  return sample;
}
