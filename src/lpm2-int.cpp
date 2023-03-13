#include <Rcpp.h>
#include "kdtree.h"
#include "uniform.h"
#include "index-list.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.lpm2_int_cpp)]]
Rcpp::IntegerVector lpm2_int_cpp(
  int n,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method
) {
  int N = x.ncol();
  double *xx = REAL(x);

  int *probability = new int[N];
  IndexList *idx = new IndexList(N);
  int *neighbours = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = n;
    idx->set(i);
  }

  while (idx->length() > 1) {
    int id1 = idx->draw();
    int len = tree->findNeighbour(neighbours, N, id1);
    int id2 = len == 1 ? neighbours[0] : neighbours[intuniform(len)];

    int p1 = probability[id1];
    int p2 = probability[id2];
    int psum = p1 + p2;

    if (psum > N) {
      if (N - p2 > intuniform((N << 1) - psum)) {
        probability[id1] = N;
        probability[id2] = psum - N;
      } else {
        probability[id1] = psum - N;
        probability[id2] = N;
      }
    } else {
      if (p2 > intuniform(psum)) {
        probability[id1] = 0;
        probability[id2] = psum;
      } else {
        probability[id1] = psum;
        probability[id2] = 0;
      }
    }

    if (probability[id1] == 0 || probability[id1] == N) {
      idx->erase(id1);
      tree->removeUnit(id1);
    }

    if (probability[id2] == 0 || probability[id2] == N) {
      idx->erase(id2);
      tree->removeUnit(id2);
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (intuniform(N) < probability[id1])
      probability[id1] = N;
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
  delete idx;
  delete tree;

  return sample;
}
