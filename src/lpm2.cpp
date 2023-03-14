#include <Rcpp.h>
#include "lpm-internal.h"

// [[Rcpp::export(.lpm2_cpp)]]
Rcpp::IntegerVector lpm2_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  double* xx = REAL(x);

  double* probabilities = new double[N];
  int* sample = new int[N];
  int sampleSize = 0;

  IndexList* idx = new IndexList(N);
  KDTree* tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    idx->set(i);
  }

  Lpm2Search* lpm2search = new Lpm2Search(tree, idx, N);

  std::function<void (int*)> unitfun = [&lpm2search](int* pair) {
    lpm2search->search(pair);
    return;
  };

  lpm_internal(
    tree,
    idx,
    probabilities,
    N,
    eps,
    sample,
    &sampleSize,
    unitfun
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);
  delete[] probabilities;
  delete[] sample;
  delete idx;
  delete tree;
  delete lpm2search;

  return svec;
}

// [[Rcpp::export(.lpm2_int_cpp)]]
Rcpp::IntegerVector lpm2_int_cpp(
  int n,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method
) {
  int N = x.ncol();
  double* xx = REAL(x);

  int* probabilities = new int[N];
  int* sample = new int[N];
  int sampleSize = 0;

  IndexList* idx = new IndexList(N);
  KDTree* tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = n;
    idx->set(i);
  }

  Lpm2Search* lpm2search = new Lpm2Search(tree, idx, N);

  std::function<void (int*)> unitfun = [&lpm2search](int* pair) {
    lpm2search->search(pair);
    return;
  };

  lpm_int_internal(
    tree,
    idx,
    probabilities,
    N,
    sample,
    &sampleSize,
    unitfun
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);
  delete[] probabilities;
  delete[] sample;
  delete idx;
  delete tree;
  delete lpm2search;

  return svec;
}
