#include <Rcpp.h>
#include "kdtree.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// [[Rcpp::export(.vsb0_cpp)]]
double vsb0_cpp(
  Rcpp::NumericVector &probs,
  Rcpp::NumericVector &ys,
  Rcpp::NumericMatrix &xs,
  int bucketSize,
  int method
) {
  int N = xs.ncol();
  double *xx = REAL(xs);
  int *neighbours = new int[N];
  double *yp = new double[N];

  KDTree *tree = new KDTree(xx, N, xs.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++)
    yp[i] = ys[i] / probs[i];

  double result = 0.0;

  for (int i = 0; i < N; i++) {
    int len = tree->findNeighbour(neighbours, N, i);

    double localMean = yp[i];

    for (int j = 0; j < len; j++)
      localMean += yp[neighbours[j]];

    localMean = yp[i] - localMean / (double)(len + 1);
    result += (double)(len + 1) / (double)(len) * (localMean * localMean);
  }

  delete[] neighbours;
  delete[] yp;
  delete tree;

  return result;
}

// [[Rcpp::export(.vsbn_cpp)]]
double vsbn_cpp(
  Rcpp::NumericVector &probs,
  Rcpp::NumericVector &ys,
  Rcpp::NumericMatrix &xs,
  int n,
  int bucketSize,
  int method
) {
  if (n < 1)
    std::range_error("n must be > 1");

  int N = xs.ncol();
  double *xx = REAL(xs);
  int *neighbours = new int[N];
  double *yp = new double[N];
  double *distances = new double[N];

  KDTree *tree = new KDTree(xx, N, xs.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++)
    yp[i] = ys[i] / probs[i];

  double result = 0.0;

  for (int i = 0; i < N; i++) {
    int len = tree->findNeighboursN(distances, neighbours, n, i);

    double localMean = yp[i];

    for (int j = 0; j < len; j++)
      localMean += yp[neighbours[j]];

    localMean = yp[i] - localMean / (double)(len + 1);
    result += (double)(len + 1) / (double)(len) * (localMean * localMean);
  }

  delete[] neighbours;
  delete[] yp;
  delete[] distances;
  delete tree;

  return result;
}


