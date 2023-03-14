#include <algorithm>
#include <Rcpp.h>
#include "scps-internal.h"

// [[Rcpp::export(.scps_coord_cpp)]]
Rcpp::IntegerVector scps_coord_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::NumericVector &random,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  double *xx = REAL(x);

  double *probabilities = new double[N];
  int *sample = new int[N];
  int sampleSize = 0;

  IndexList *idx = new IndexList(N);
  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  std::function<double (const int)> randfun = [&random](int i) {
    return random[i];
  };

  int id = 0;
  std::function<int ()> unitfun = [&idx, &id]() {
    while(!idx->exists(id))
      id += 1;

    int tid = id;
    id += 1;

    return tid;
  };

  scps_internal(
    tree,
    idx,
    probabilities,
    N,
    eps,
    sample,
    &sampleSize,
    randfun,
    unitfun
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete[] sample;
  delete idx;
  delete tree;

  return svec;
}
