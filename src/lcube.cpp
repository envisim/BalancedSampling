#include <algorithm>
#include <Rcpp.h>
#include "cube-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))
#define mati(r, c, p) ((r) * (p) + (c))

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
  int p = xbal.ncol();
  int sampleSize = 0;

  if (N != xspread.ncol())
    std::range_error("xbal != xspread");

  double *amat = new double[N * p];
  double *probabilities = new double[N];
  IndexList *idx = new IndexList(N);
  int *sample = new int[N];

  double *xxspread = REAL(xspread);
  KDTree *tree = new KDTree(xxspread, N, xspread.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  for (int i = 0; i < N; i++) {
    if (pclose(prob[i], eps)) {
      if (pbig(prob[i], eps)) {
        sample[sampleSize] = i;
        sampleSize += 1;
      }

      idx->erase(i);
      tree->removeUnit(i);

      continue;
    }

    for (int k = 0; k < p; k++) {
      amat[mati(k, i, N)] = xbal(i, k) / prob[i];
    }
  }

  cubeRunPhases(
    amat,
    probabilities,
    idx,
    tree,
    sample,
    &sampleSize,
    eps,
    p,
    N
  );

  std::sort(sample, sample + sampleSize);
  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] amat;
  delete[] probabilities;
  delete idx;
  delete[] sample;
  delete tree;

  return svec;
}
