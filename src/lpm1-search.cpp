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

void traverse(
  int *idx1,
  int *idx2,
  KDTree* tree,
  std::unordered_set<int> &idx,
  int *neighbours1,
  int *neighbours2,
  int *history,
  int *histn,
  int N
) {
  *idx1 = history[*histn - 1];
  int len1 = tree->findNeighbour(neighbours1, N, *idx1);
  int len1_2 = len1;

  for (int i = 0; i < len1;) {
    int len2 = tree->findNeighbour(neighbours2, N, neighbours1[i]);
    int found = 0;
    for (int j = 0; j < len2; j++) {
      if (neighbours2[j] == *idx1) {
        found = 1;
        break;
      }
    }

    if (found == 0) {
      len1 -= 1;
      if (i != len1) {
        int temp = neighbours1[i];
        neighbours1[i] = neighbours1[len1];
        neighbours1[len1] = temp;
      }
    } else {
      i += 1;
    }
  }

  if (len1 > 0) {
    *idx2 = len1 == 1 ? neighbours1[0] : neighbours1[intuniform(len1)];
    return;
  }

  if (*histn == N) {
    *histn = 0;
  }

  history[*histn] = len1_2 == 1 ? neighbours1[0] : neighbours1[intuniform(len1_2)];
  *histn += 1;
  traverse(idx1, idx2, tree, idx, neighbours1, neighbours2, history, histn, N);
}

void search(
  int *idx1,
  int *idx2,
  KDTree* tree,
  std::unordered_set<int> &idx,
  int *neighbours1,
  int *neighbours2,
  int *history,
  int *histn,
  int N
) {
  while (*histn > 0) {
    if (idx.count(history[*histn - 1])) {
      break;
    }

    *histn -= 1;
  }

  if (*histn == 0) {
    int u1 = intuniform(N);
    std::unordered_set<int>::iterator it1 = idx.find(u1);
    bool exists = it1 != idx.end();
    if (exists) {
      history[0] = u1;
    } else {
      idx.insert(u1);
      it1 = idx.find(u1);
      it1++;
      history[0] = it1 == idx.end() ? *idx.begin() : *it1;
      idx.erase(u1);
    }

    *histn = 1;
  }

  traverse(idx1, idx2, tree, idx, neighbours1, neighbours2, history, histn, N);
}

// [[Rcpp::export(.lpm1_search_cpp)]]
Rcpp::IntegerVector lpm1_search_cpp(
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
  int *neighbours1 = new int[N];
  int *neighbours2 = new int[N];
  int *history = new int[N];
  int histn = 0;


  KDTree *tree = new KDTree(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probability[i] = prob[i];
    idx.insert(i);
  }

  while (unresolvedObjects > 1) {
    int idx1;
    int idx2;
    search(
      &idx1,
      &idx2,
      tree,
      idx,
      neighbours1,
      neighbours2,
      history,
      &histn,
      N
    );

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
      histn -= 1;
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
      neighbours1[j] = i + 1;
      j += 1;
    }
  }

  Rcpp::IntegerVector sample(neighbours1, neighbours1 + j);
  delete[] neighbours1;
  delete[] neighbours2;
  delete[] probability;
  delete[] history;
  delete tree;

  return sample;
}
