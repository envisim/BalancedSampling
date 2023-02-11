#include <Rcpp.h>
#include "kdtree.h"

//**********************************************
// Author: Wilmer Prentius
// Last edit: 2023-02-10
// Licence: GPL (>=2)
//**********************************************

// struct Object {
//   int probability;
//   unsigned int index;
//   Object() : probability(0), index(0) {}
// };

int randn(double u, int n) {
  return (int)((double)n * u);
}

// [[Rcpp::export]]
Rcpp::IntegerVector lpm2inttree(int n, Rcpp::NumericMatrix &x) {
  int N = x.ncol();
  int unresolvedObjects = N;
  double *xx = REAL(x);

  int *probability = new int[N];
  int *idx = new int[N];
  int *neighbours = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), 40);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = n;
    idx[i] = i;
  }

  Rcpp::NumericVector rand1 = Rcpp::runif(N, 0.0, 1.0);
  Rcpp::NumericVector rand2 = Rcpp::runif(N, 0.0, 1.0);

  while (unresolvedObjects > 1) {
    int u1 = randn(rand2[unresolvedObjects-1], unresolvedObjects);
    int idx1 = idx[u1];
    int len = tree->findNeighbour(neighbours, N, idx1);
    int idx2 = len == 1 ? neighbours[0] : randn(R::runif(0.0, 1.0), len);

    int p1 = probability[idx1];
    int p2 = probability[idx2];
    int psum = p1 + p2;
    double u = rand1[unresolvedObjects-1];

    if (psum > N) {
      if (N - p2 > randn(u, (N << 1) - psum)) {
        probability[idx1] = N;
        probability[idx2] = psum - N;
      } else {
        probability[idx1] = psum - N;
        probability[idx2] = N;
      }
    } else {
      if (p2 > randn(u, psum)) {
        probability[idx1] = 0;
        probability[idx2] = psum;
      } else {
        probability[idx1] = psum;
        probability[idx2] = 0;
      }
    }

    if (probability[idx2] == 0 || probability[idx2] == N) {
      int u2;
      for (int i = 0; i < unresolvedObjects; i++) {
        if (idx[i] == idx2) {
          u2 = i;
          break;
        }
      }

      unresolvedObjects -= 1;

      int temp = idx[unresolvedObjects];
      idx[unresolvedObjects] = idx2;
      idx[u2] = temp;
      tree->removeUnit(idx2);

      if (unresolvedObjects == u1)
        u1 = u2;
    }

    if (probability[idx1] == 0 || probability[idx1] == N) {
      unresolvedObjects -= 1;
      int temp = idx[unresolvedObjects];
      idx[unresolvedObjects] = idx1;
      idx[u1] = temp;
      tree->removeUnit(idx1);
    }
  }

  delete[] idx;
  delete[] neighbours;
  delete[] tree->top->units;
  delete tree->top;
  delete tree;
  Rcpp::IntegerVector s(n);

  for (int i = 0, j = 0; i < N && j < n; i++) {
    if (probability[i] == N) {
      s[j] = i + 1;
      j += 1;
    }
  }

  delete[] probability;

  return s;
}
