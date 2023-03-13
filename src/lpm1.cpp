#include <Rcpp.h>
#include "kdtree.h"
#include "uniform.h"
#include "index-list.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))

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

  double *probability = new double[N];
  IndexList *idx = new IndexList(N);
  int *neighbours1 = new int[N];
  int *neighbours2 = new int[N];

  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = prob[i];
    idx->set(i);
  }

  while (idx->length() > 1) {
    int id1 = idx->draw();
    int len1 = tree->findNeighbour(neighbours1, N, id1);

    for (int i = 0; i < len1;) {
      int len2 = tree->findNeighbour(neighbours2, N, neighbours1[i]);
      int found = 0;
      for (int j = 0; j < len2; j++) {
        if (neighbours2[j] == id1) {
          found = 1;
          break;
        }
      }

      if (found == 0) {
        len1 -= 1;
        neighbours1[i] = neighbours1[len1];
      } else {
        i += 1;
      }
    }

    if (len1 < 1)
      continue;

    int id2 = len1 == 1 ? neighbours1[0] : neighbours1[intuniform(len1)];

    double p1 = probability[id1];
    double p2 = probability[id2];
    double psum = p1 + p2;

    if (psum > 1.0) {
      if (1.0 - p2 > stduniform() * (2.0 - psum)) {
        probability[id1] = 1.0;
        probability[id2] = psum - 1.0;
      } else {
        probability[id1] = psum - 1.0;
        probability[id2] = 1.0;
      }
    } else {
      if (p2 > stduniform() * psum) {
        probability[id1] = 0.0;
        probability[id2] = psum;
      } else {
        probability[id1] = psum;
        probability[id2] = 0.0;
      }
    }

    if (pclose(probability[id1], eps)) {
      idx->erase(id1);
      tree->removeUnit(id1);
    }

    if (pclose(probability[id2], eps)) {
      idx->erase(id2);
      tree->removeUnit(id2);
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (stduniform() < probability[id1])
      probability[id1] = 1.0;
  }

  int j = 0;
  for (int i = 0; i < N; i++) {
    if (probability[i] >= 1.0 - eps) {
      neighbours1[j] = i + 1;
      j += 1;
    }
  }

  Rcpp::IntegerVector sample(neighbours1, neighbours1 + j);
  delete[] probability;
  delete[] neighbours1;
  delete[] neighbours2;
  delete idx;
  delete tree;

  return sample;
}

