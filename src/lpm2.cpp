#include "lpm2-internal.h"
#include <Rcpp.h>

// [[Rcpp::export(.lpm2_cpp)]]
Rcpp::IntegerVector lpm2_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  double *xx = REAL(x);

  double *probabilities = new double[N];
  int *neighbours = new int[N];
  int *sample = new int[N];
  int sampleSize = 0;

  IndexList *idx = new IndexList(N);
  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    idx->set(i);
  }

  lpm2_internal(
    tree,
    idx,
    probabilities,
    neighbours,
    N,
    eps,
    sample,
    &sampleSize
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);
  delete[] probabilities;
  delete[] neighbours;
  delete[] sample;
  delete idx;
  delete tree;

  return svec;
}
