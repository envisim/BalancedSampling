#include <Rcpp.h>
#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// struct Lpm1Search : Lpm2Search {
//   int* tneighbours;
//   Lpm1Search(KDTree* t_tree, IndexList* t_idx, const int t_N) {
//     Lpm2Search
//   }
// };

// [[Rcpp::export(.lpm1_cpp)]]
Rcpp::IntegerVector lpm1_cpp(
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
  int *tneighbours = new int[N];
  int *sample = new int[N];
  int sampleSize = 0;

  IndexList *idx = new IndexList(N);
  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    idx->set(i);
  }

  std::function<void (int*)> unitfun = [&tree, &idx, &neighbours, &tneighbours, &N](int* pair) {
    while (true) {
      pair[0] = idx->draw();
      int len = tree->findNeighbour(neighbours, N, pair[0]);

      for (int i = 0; i < len;) {
        int tlen = tree->findNeighbour(tneighbours, N, neighbours[i]);
        bool found = false;

        for (int j = 0; j < tlen; j++) {
          if (tneighbours[j] == pair[0]) {
            found = true;
            break;
          }
        }

        if (found) {
          i += 1;
        } else {
          len -= 1;
          neighbours[i] = neighbours[len];
        }
      }

      if (len > 0) {
        pair[1] = neighbours[intuniform(len)];
        return;
      }
    }
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
  delete[] neighbours;
  delete[] tneighbours;
  delete[] sample;
  delete idx;
  delete tree;

  return svec;
}
