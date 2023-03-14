#include <Rcpp.h>
#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

struct Lpm1Search {
  KDTree* tree;
  IndexList* idx;
  int* neighbours;
  int* tneighbours;
  int* history;
  int histn;
  int N;
  Lpm1Search(KDTree* t_tree, IndexList* t_idx, const int t_N) {
    tree = t_tree;
    idx = t_idx;
    neighbours = new int[t_N];
    tneighbours = new int[t_N];
    history = new int[t_N];
    histn = 0;
    N = t_N;
  };
  ~Lpm1Search() {
    delete[] neighbours;
    delete[] tneighbours;
    delete[] history;
  };
  void traverse(int* pair) {
    // Set the first unit to the last in history
    pair[0] = history[histn - 1];
    // Find this units nearest neighbours
    int len = tree->findNeighbour(neighbours, N, pair[0]);
    int len_copy = len;

    // Go through all nearest neighbours
    for (int i = 0; i < len;) {
      // Find the neighbours nearest neighbours
      int tlen = tree->findNeighbour(tneighbours, N, neighbours[i]);
      bool found = false;

      // Check if any of these are the history-unit
      for (int j = 0; j < tlen; j++) {
        if (tneighbours[j] == pair[0]) {
          found = true;
          break;
        }
      }

      // If the history-unit exists among the nearest neighbours, we continue
      // to see if any other of the history-units neighbours also are mutual.
      // Otherwise, the history-unit is not among the nearest neighbours,
      // we swap places and continue the search.
      if (found) {
        i += 1;
      } else {
        len -= 1;
        if (i != len) {
          int temp = neighbours[i];
          neighbours[i] = neighbours[len];
          neighbours[len] = temp;
        }
      }
    }

    // If we found one or more mutual neighbours, we select one at random
    if (len > 0) {
      pair[1] = neighbours[intuniform(len)];
      return;
    }

    // If we come here, no mutual neighbours exist

    // We might need to clear the history if the search has been going on for
    // too long. This can probably? happen if there is a long history, and
    // updates has affected previous units.
    if (histn == N)
      histn = 0;

    // We select a unit at random to become the next history unit, and traverse
    // one step further.
    history[histn] = neighbours[intuniform(len_copy)];
    histn += 1;
    traverse(pair);
    return;
  }
  void search(int* pair) {
    // Go back in the history and remove units that does not exist
    while (histn > 0) {
      if (idx->exists(history[histn - 1])) {
        break;
      }

      histn -= 1;
    }

    // If there is no history, we draw a unit at random
    if (histn == 0) {
      history[0] = idx->draw();
      histn = 1;
    }

    // Traverse the history to find a unit with a mutual nearest neighbour
    traverse(pair);
    return;
  };
};


// [[Rcpp::export(.lpm1_search_cpp)]]
Rcpp::IntegerVector lpm1_search_cpp(
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

  Lpm1Search* lpm1search = new Lpm1Search(tree, idx, N);

  std::function<void (int*)> unitfun = [&lpm1search](int* pair) {
    lpm1search->search(pair);
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
  delete lpm1search;

  return svec;
}
