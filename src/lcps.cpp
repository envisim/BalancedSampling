#include <algorithm>
#include <Rcpp.h>
#include "uniform.h"
#include "scps-internal.h"

// [[Rcpp::export(.lcps_cpp)]]
Rcpp::IntegerVector lcps_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  double *xx = REAL(x);

  double *probabilities = new double[N];
  int *sample = new int[N];
  int sampleSize = 0;

  IndexList *idx = new IndexList(N);
  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  // Prepare arrays needed by the findNeighbours algorithm
  double *idxarr = new double[N];
  double *dists = new double[N];
  double *weights = new double[N];
  int *neighbours = new int[N];

  std::function<double (const int)> randfun = [](int i) {
    return stduniform();
  };

  // The function to choose the next deciding unit
  std::function<int ()> unitfun = [&idx, &tree, &probabilities, &idxarr, &dists, &weights, &neighbours]() {
    // Take care of edge cases
    if (idx->length() <= 1) {
      if (idx->length() == 1)
        return idx->get(0);
      if (idx->length() < 1)
        std::range_error("trying to find index in empty list");
    }

    int idxarrlen = 0;
    double idxmindist = DBL_MAX;

    // Loop through all remaining units.
    // Put the smallest distances in idxarr
    for (int i = 0; i < idx->length(); i++) {
      int len = tree->findNeighbours(probabilities, weights, dists, neighbours, idx->get(i));
      double idxdist = dists[neighbours[len - 1]];

      if (idxdist < idxmindist) {
        idxarr[0] = i;
        idxmindist = idxdist;
        idxarrlen = 1;
      } else if (idxdist == idxmindist) {
        idxarr[idxarrlen] = i;
        idxarrlen += 1;
      }
    }

    // Choose randomly from the units in idxarr
    int idxi = intuniform(idxarrlen);
    return idx->get(idxi);
  };

  scps_internal(
    tree,
    idx,
    probabilities,
    N,
    eps,
    sample,
    &sampleSize,
    randfun,
    unitfun
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete[] sample;
  delete idx;
  delete tree;
  delete[] idxarr;
  delete[] dists;
  delete[] weights;
  delete[] neighbours;

  return svec;
}
