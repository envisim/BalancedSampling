#include <Rcpp.h>
#include "lpm-internal.h"

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
  int p = x.nrow();
  double* xx = REAL(x);

  if (prob.length() != N)
    std::invalid_argument("prob an x does not match");

  // int *sample = new int[N];
  // int sampleSize = 0;

  KDTree* orgTree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  orgTree->init();
  KDTree* tree = orgTree->copy();

  IndexList* orgIdx = new IndexList(N);
  double* probabilities = new double[N];

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    orgIdx->set(i);
  }

  Lpm lpm(probabilities, tree, orgIdx, N, LpmMethod::lpm2, eps);
  lpm.run();
  int sampleSize = lpm.sampleSize;

  delete tree;
  orgIdx->reset();

  for (int i = 0, j = 0; i < N; i++) {
    if (j < sampleSize && i == lpm.sample[j] - 1) {
      j += 1;
      continue;
    }

    orgTree->removeUnit(i);
    orgIdx->erase(i);
  }

  // orgTree->prune();

  Rcpp::IntegerMatrix smat(sampleSize, 2);
  int* psmat = INTEGER(smat);

  for (int i = 0; i < sampleSize; i++)
    psmat[i] = lpm.sample[i];

  int remainingSize = sampleSize;

  for (int i = 0; i < sn - 1; i++) {
    double subprob = (double)sizes[i] / (double)remainingSize;
    KDTree* tree = orgTree->copy();
    IndexList* idx = orgIdx->copyLen();

    for (int j = 0; j < idx->length(); j++) {
      probabilities[idx->get(j)] = subprob;
    }

    lpm.probabilities = probabilities;
    lpm.tree = tree;
    lpm.idx = idx;
    lpm.sampleSize = 0;
    lpm.run();

    for (int j = 0, k = 0; j < sampleSize && k < lpm.sampleSize; j++) {
      if (psmat[j] == lpm.sample[k]) {
        orgTree->removeUnit(lpm.sample[k] - 1);
        orgIdx->erase(lpm.sample[k] - 1);
        psmat[j + sampleSize] = i + 1;
        k += 1;
      }
    }

    remainingSize -= lpm.sampleSize;
    delete tree;
    delete idx;
  }

  // Set the remaining units to the remaining group
  for (int i = 0; i < sampleSize; i++) {
    if (psmat[sampleSize + i] == 0)
      psmat[sampleSize + i] = sn;
  }

  delete[] probabilities;
  delete orgIdx;
  delete orgTree;

  return smat;
}

