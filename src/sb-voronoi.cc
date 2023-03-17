#include <algorithm>
#include <Rcpp.h>
#include "kdtree.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

//[[Rcpp::export(.sb_voronoi_cpp)]]
double sb_voronoi_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &sample,
  int bucketSize,
  int method
) {
  int N = x.ncol();
  int n = sample.length();
  int p = x.nrow();
  double *xx = REAL(x);
  double *xs = new double[n*p];
  double *incl = new double[n];
  int *neighbours = new int[n];

  for (int i = 0; i < n; i++) {
    if (sample[i] < 1 || sample[i] > N)
      Rcpp::stop("'sample' must contain unit indices");

    std::copy_n(xx + (sample[i] - 1) * p, p, xs + i * p);
    incl[i] = 0.0;
  }

  KDTree *tree = new KDTree(xs, n, p, bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    int len = tree->findClosest(neighbours, n, xx + i * p);

    double ppart = len == 1 ? prob[i] : prob[i] / (double)len;
    for (int j = 0; j < len; j++)
      incl[neighbours[j]] += ppart;
  }

  double result = 0.0;
  for (int i = 0; i < n; i++) {
    double temp = incl[i] - 1.0;
    result += temp * temp;
  }

  delete[] xs;
  delete[] incl;
  delete[] neighbours;
  delete tree;

  return result / (double)n;
}
