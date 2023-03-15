#include <Rcpp.h>
#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

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

  Lpm lpm(probabilities, nullptr, idx, N, rpm, eps);
  lpm.run();

  Rcpp::IntegerVector sample(lpm.sample, lpm.sample + lpm.sampleSize);

  delete idx;
  delete[] probabilities;

  return sample;
}
