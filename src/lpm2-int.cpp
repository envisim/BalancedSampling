#include <unordered_set>
#include <Rcpp.h>
#include "kdtree.h"
#include "uniform.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define intuniform(N) ((int)((double)N * stduniform()))

// [[Rcpp::export(.lpm2_int_cpp)]]
Rcpp::IntegerVector lpm2_int_cpp(
  int n,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method
) {
  int N = x.ncol();
  int unresolvedObjects = N;
  double *xx = REAL(x);

  int *probability = new int[N];
  std::unordered_set<int> idx(N);
  int *neighbours = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = n;
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
    // int idx2 = len == 1 ? neighbours[0] : neighbours[(int)((double)len * R::runif(0.0, 1.0))];
    int idx2 = len == 1 ? neighbours[0] : neighbours[intuniform(len)];

    int p1 = probability[idx1];
    int p2 = probability[idx2];
    int psum = p1 + p2;

    if (psum > N) {
      if (N - p2 > intuniform((N << 1) - psum)) {
        probability[idx1] = N;
        probability[idx2] = psum - N;
      } else {
        probability[idx1] = psum - N;
        probability[idx2] = N;
      }
    } else {
      if (p2 > intuniform(psum)) {
        probability[idx1] = 0;
        probability[idx2] = psum;
      } else {
        probability[idx1] = psum;
        probability[idx2] = 0;
      }
    }

    if (probability[idx2] == 0 || probability[idx2] == N) {
      unresolvedObjects -= 1;
      idx.erase(idx2);
      tree->removeUnit(idx2);
    }

    if (probability[idx1] == 0 || probability[idx1] == N) {
      unresolvedObjects -= 1;
      idx.erase(idx1);
      tree->removeUnit(idx1);
    }
  }

  if (unresolvedObjects == 1) {
    int idx1 = *idx.begin();
    if (intuniform(N) < probability[idx1])
      probability[idx1] = 1.0;
  }

  Rcpp::IntegerVector sample(n);
  for (int i = 0, j = 0; i < N && j < n; i++) {
    if (probability[i] == N) {
      sample[j] = i + 1;
      j += 1;
    }
  }

  delete[] probability;
  delete[] neighbours;
  delete tree;

  return sample;
}
