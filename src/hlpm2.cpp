#include <Rcpp.h>
#include "lpm-internal.h"

struct Hlpm2Search {
  KDTree* tree;
  IndexList* idx;
  int* neighbours;
  int N;
  Hlpm2Search(KDTree* t_tree, IndexList* t_idx, const int t_N) {
    tree = t_tree;
    idx = t_idx;
    neighbours = new int[t_N];
    N = t_N;
  };
  ~Hlpm2Search() {
    delete[] neighbours;
  };
  void search(int* pair) {
    pair[0] = idx->draw();
    int len = tree->findNeighbour(neighbours, N, pair[0]);
    pair[1] = neighbours[intuniform(len)];
    return;
  };
};

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
  // int *neighbours = new int[N];

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    orgIdx->set(i);
  }

  Hlpm2Search* hlpm2search = new Hlpm2Search(tree, orgIdx, N);

  std::function<void (int*)> unitfun = [&hlpm2search](int *pair) {
    hlpm2search->search(pair);
    return;
  };

  lpm_internal(
    tree,
    orgIdx,
    probabilities,
    N,
    eps,
    sample,
    &sampleSize,
    unitfun
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
    hlpm2search->tree = tree;
    hlpm2search->idx = idx;

    lpm_internal(
      tree,
      idx,
      probabilities,
      N,
      eps,
      sample,
      &rsize,
      unitfun
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
  // delete[] neighbours;
  delete[] sample;
  delete hlpm2search;
  delete orgIdx;
  delete orgTree;

  return smat;
}

