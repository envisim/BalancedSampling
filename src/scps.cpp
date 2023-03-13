#include <algorithm>
#include <Rcpp.h>
#include "uniform.h"
#include "scps-internal.h"

// [[Rcpp::export(.scps_cpp)]]
Rcpp::IntegerVector scps_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  double *xx = REAL(x);

  double *probabilities = new double[N];
  int *sample = new int[N];
  int sampleSize = 0;

  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
  }

  std::function<double (const int)> randfun = [](int i) {
    return stduniform();
  };

  scps_internal(
    xx,
    N,
    probabilities,
    tree,
    eps,
    randfun,
    sample,
    &sampleSize
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete[] sample;
  delete tree;

  return svec;
}

