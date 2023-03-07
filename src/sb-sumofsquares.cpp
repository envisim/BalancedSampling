#include <Rcpp.h>
#include "kdtree.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

//[[Rcpp::export(.sb_sumofsquares_cpp)]]
double sb_sumofsquares_cpp(
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

  for (int i = 0; i < n; i++) {
    if (sample[i] < 1 || sample[i] > N)
      Rcpp::stop("'sample' must contain unit indices");

    std::copy_n(xx + (sample[i] - 1) * p, p, xs + i * p);
  }

  KDTree *tree = new KDTree(xs, n, p, bucketSize, method);
  tree->init();

  double result = 0.0;
  double total = 0.0;
  double *means = new double[p];

  for (int i = 0; i < N; i++) {
    double *unit = xx + i * p;
    double dist = tree->findSmallestDistanceToPoint(unit);
    result += dist;

    for (int k = 0; k < p; k++) {
      means[k] += *(unit + k);
      total += *(unit + k) * *(unit + k);
    }
  }

  for (int k = 0; k < p; k++)
    total -= means[k] * means[k] / (double)N;

  result *= (double)n / total;

  delete[] xs;
  delete[] means;
  delete tree;

  return result;
}
