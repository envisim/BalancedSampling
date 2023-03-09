#include "lpm2-internal.h"
#include <Rcpp.h>

// [[Rcpp::export(.hlpm2_cpp)]]
Rcpp::IntegerMatrix hlpm2_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &sizes,
  int bucketSize,
  int method,
  double eps
) {
  int sn = sizes.length();
  int N = x.ncol();
  double *xx = REAL(x);

  int *sample = new int[N];
  int sampleSize = 0;

  KDTree *orgTree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  orgTree->init();
  KDTree *tree = orgTree->copy();

  IndexList *orgIdx = new IndexList(N);
  double *probabilities = new double[N];
  int *neighbours = new int[N];

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    orgIdx->set(i);
  }

  lpm2_internal(
    tree,
    orgIdx,
    probabilities,
    neighbours,
    N,
    eps,
    sample,
    &sampleSize
  );

  delete tree;

  orgIdx->reset();

  for (int i = 0, j = 0; i < N; i++) {
    if (j < sampleSize && i == sample[j] - 1) {
      j += 1;
      continue;
    }

    orgTree->removeUnit(i);
    orgIdx->erase(i);
  }

  // orgTree->prune();

  Rcpp::IntegerMatrix smat(sampleSize, 2);
  int *psmat = INTEGER(smat);

  for (int i = 0; i < sampleSize; i++)
    psmat[i] = sample[i];

  int remainingSize = sampleSize;

  for (int i = 0; i < sn - 1; i++) {
    double subprob = (double)sizes[i] / (double)remainingSize;
    KDTree *tree = orgTree->copy();
    IndexList* idx = orgIdx->copyLen();

    for (int j = 0; j < idx->length(); j++) {
      probabilities[idx->get(j)] = subprob;
    }

    int rsize = 0;

    lpm2_internal(
      tree,
      idx,
      probabilities,
      neighbours,
      N,
      eps,
      sample,
      &rsize
    );

    for (int j = 0, k = 0; j < sampleSize && k < rsize; j++) {
      if (psmat[j] == sample[k]) {
        orgTree->removeUnit(sample[k] - 1);
        orgIdx->erase(sample[k] - 1);
        psmat[j + sampleSize] = i + 1;
        k += 1;
      }
    }

    remainingSize -= rsize;
    delete tree;
    delete idx;
  }

  for (int i = 0; i < sampleSize; i++) {
    if (psmat[sampleSize + i] == 0)
      psmat[sampleSize + i] = sn;
  }

  delete[] probabilities;
  delete[] neighbours;
  delete[] sample;
  delete orgIdx;
  delete orgTree;

  return smat;
}

