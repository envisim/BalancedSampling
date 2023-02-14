#include <unordered_set>
#include <Rcpp.h>
#include "kdtree.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.lpm1_cpp)]]
Rcpp::IntegerVector lpm1_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method
) {
  int N = x.ncol();
  int unresolvedObjects = N;
  double *xx = REAL(x);

  double nprob = 0.0;
  double *probability = new double[N];
  std::unordered_set<int> idx(N);
  int *neighbours = new int[N];
  int *neighbours2 = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = prob[i];
    nprob += probability[i];
    idx.insert(i);
  }

  Rcpp::NumericVector rand1 = Rcpp::runif(N, 0.0, 1.0);

  while (unresolvedObjects > 1) {
    int u1 = (int)(R::runif(0.0, 1.0) * (double)N);
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

    for (int i = 0; i < len;) {
      int len2 = tree->findNeighbour(neighbours2, N, neighbours[i]);
      int found = 0;
      for (int j = 0; j < len2; j++) {
        if (neighbours2[j] == idx1) {
          found = 1;
          break;
        }
      }

      if (found == 0) {
        neighbours[i] = neighbours[len-1];
        len -= 1;
      } else {
        i += 1;
      }
    }

    if (len < 1)
      continue;

    int idx2 = len == 1 ? neighbours[0] : (int)(R::runif(0.0, 1.0) * (double)len);

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

