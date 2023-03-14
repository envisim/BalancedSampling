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

  IndexList *idx = new IndexList(N);
  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  std::function<double (const int)> randfun = [](int i) {
    return stduniform();
  };

  std::function<int ()> unitfun = [&idx]() {
    return idx->draw();
  };

  scps_internal(
    N,
    probabilities,
    tree,
    eps,
    randfun,
    unitfun,
    sample,
    &sampleSize,
    idx
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete[] sample;
  delete idx;
  delete tree;

  return svec;
}

