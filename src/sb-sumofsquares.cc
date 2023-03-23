#include <algorithm>
#include <stddef.h>

#include <Rcpp.h>

#include "KDStoreClass.h"
#include "KDTreeClass.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

//[[Rcpp::export(.sb_sumofsquares_cpp)]]
double sb_sumofsquares_cpp(
  Rcpp::NumericMatrix &x,
  Rcpp::IntegerVector &sample,
  size_t treeBucketSize,
  size_t treeMethod
) {
  size_t N = x.ncol();
  size_t p = x.nrow();
  size_t n = sample.length();
  double* xx = REAL(x);
  double* xs = new double[n * p];

  for (size_t i = 0; i < n; i++) {
    if (sample[i] < 1 || sample[i] > N)
      throw std::range_error("'sample' must contain unit indices");

    std::copy_n(xx + (sample[i] - 1) * p, p, xs + i * p);
  }

  KDTree tree(xs, n, p, treeBucketSize, IntToKDTreeSplitMethod(treeMethod));
  KDStore store(n, 1);

  double result = 0.0;
  double total = 0.0;
  double* means = new double[p];

  for (size_t i = 0; i < N; i++) {
    double* unit = xx + i * p;
    tree.FindNeighbours(&store, unit);
    double dist = store.MinimumDistance();
    result += dist;

    for (size_t k = 0; k < p; k++) {
      means[k] += *(unit + k);
      total += *(unit + k) * *(unit + k);
    }
  }

  for (size_t k = 0; k < p; k++)
    total -= means[k] * means[k] / (double)N;

  result *= (double)n / total;

  delete[] xs;
  delete[] means;

  return result;
}
