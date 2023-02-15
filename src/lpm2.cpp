#include <unordered_set>
#include <Rcpp.h>
#include "kdtree.h"
#include "uniform.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define intuniform(N) ((int)((double)N * stduniform()))
#define pclose(p, eps) (p <= eps || p >= 1.0 - eps)

// [[Rcpp::export(.lpm2_cpp)]]
Rcpp::IntegerVector lpm2_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  int unresolvedObjects = N;
  double *xx = REAL(x);

  double *probability = new double[N];
  std::unordered_set<int> idx(N);
  int *neighbours = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = prob[i];
    idx.insert(i);
  }

  while (unresolvedObjects > 1) {
    int u1 = intuniform(N);
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
    int idx2 = len == 1 ? neighbours[0] : neighbours[intuniform(len)];

    double p1 = probability[idx1];
    double p2 = probability[idx2];
    double psum = p1 + p2;

    if (psum > 1.0) {
      if (1.0 - p2 > stduniform() * (2.0 - psum)) {
        probability[idx1] = 1.0;
        probability[idx2] = psum - 1.0;
      } else {
        probability[idx1] = psum - 1.0;
        probability[idx2] = 1.0;
      }
    } else {
      if (p2 > stduniform() * psum) {
        probability[idx1] = 0.0;
        probability[idx2] = psum;
      } else {
        probability[idx1] = psum;
        probability[idx2] = 0.0;
      }
    }

    if (pclose(probability[idx2], eps)) {
      unresolvedObjects -= 1;
      idx.erase(idx2);
      tree->removeUnit(idx2);
    }

    if (pclose(probability[idx1], eps)) {
      unresolvedObjects -= 1;
      idx.erase(idx1);
      tree->removeUnit(idx1);
    }
  }

  if (unresolvedObjects == 1) {
    int idx1 = *idx.begin();
    if (stduniform() < probability[idx1])
      probability[idx1] = 1.0;
  }

  int j = 0;
  for (int i = 0; i < N; i++) {
    if (probability[i] == 1.0) {
      neighbours[j] = i + 1;
      j += 1;
    }
  }

  Rcpp::IntegerVector sample(neighbours, neighbours + j);
  delete[] probability;
  delete[] neighbours;
  delete tree;

  return sample;
}
