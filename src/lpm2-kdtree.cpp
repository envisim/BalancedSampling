#include <Rcpp.h>
#include <unordered_set>
#include "kdtree.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export]]
Rcpp::IntegerVector lpm2tree(Rcpp::NumericVector &prob, Rcpp::NumericMatrix &x) {
  int N = x.ncol();
  int unresolvedObjects = N;
  double *xx = REAL(x);

  double nprob = 0.0;
  double *probability = new double[N];
  // int *idx = new int[N];
  std::unordered_set<int> idx(N);
  int *neighbours = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), 40, 2);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = prob[i];
    nprob += probability[i];
    // idx[i] = i;
    idx.insert(i);
  }

  Rcpp::NumericVector rand1 = Rcpp::runif(N, 0.0, 1.0);
  Rcpp::NumericVector rand2 = Rcpp::runif(N, 0.0, 1.0);

  while (unresolvedObjects > 1) {
    int u1 = rand2[unresolvedObjects-1] * (double)N;
    std::unordered_set<int>::iterator it1 = idx.find(u1);
    bool exists = it1 != idx.end();
    int idx1;
    if (exists) {
      idx1 = u1;
    } else {
      idx.insert(u1);
      it1 = idx.find(u1);
      it1++;
      idx1 = it1 == idx.end() ? *idx.begin() : *it1;
      idx.erase(u1);
    }

    int len = tree->findNeighbour(neighbours, N, idx1);
    int idx2 = len == 1 ? neighbours[0] : (int)((double)len * R::runif(0.0, 1.0));

    double p1 = probability[idx1];
    double p2 = probability[idx2];
    double psum = p1 + p2;
    double u = rand1[unresolvedObjects-1];

    if (psum > 1.0) {
      if (1.0 - p2 > u * (2.0 - psum)) {
        probability[idx1] = 1.0;
        probability[idx2] = psum - 1.0;
      } else {
        probability[idx1] = psum - 1.0;
        probability[idx2] = 1.0;
      }
    } else {
      if (p2 > u * psum) {
        probability[idx1] = 0.0;
        probability[idx2] = psum;
      } else {
        probability[idx1] = psum;
        probability[idx2] = 0.0;
      }
    }

    if (probability[idx2] == 0.0 || probability[idx2] == 1.0) {
      unresolvedObjects -= 1;
      idx.erase(idx2);
      tree->removeUnit(idx2);
    }

    if (probability[idx1] == 0.0 || probability[idx1] == 1.0) {
      unresolvedObjects -= 1;
      idx.erase(idx1);
      tree->removeUnit(idx1);
    }
  }

  if (unresolvedObjects == 1) {
    int idx1 = *idx.begin();
    if (rand1[0] < probability[idx1])
      probability[idx1] = 1.0;
  }

  delete[] neighbours;
  delete tree;
  int n = (int)nprob + 1;
  int j = 0;
  int *s = new int[n];

  for (int i = 0; i < N && j < n; i++) {
    if (probability[i] == 1.0) {
      s[j] = i + 1;
      j += 1;
    }
  }

  delete[] probability;

  Rcpp::IntegerVector s2(s, s + j);
  delete[] s;

  return s2;
}
