#include <algorithm>
#include <functional>
#include <Rcpp.h>
#include "kdtree-cps.h"
#include "uniform.h"
#include "index-list.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

// #define intuniform(N) ((int)((double)N * stduniform()))
// #define pclose(p, eps) (p <= eps || p >= 1.0 - eps)

void decide(
  IndexList *idx,
  KDTreeCps *tree,
  double *probs,
  int *sampleSize,
  int uid,
  double eps
) {
  if (probs[uid] <= eps) {
    idx->erase(uid);
    tree->removeUnit(uid);
  } else if (probs[uid] >= 1.0 - eps) {
    idx->erase(uid);
    tree->removeUnit(uid);
    *sampleSize += 1;
  }
}

void scps_internal(
  const double *xx,
  const int N,
  double *probabilities,
  KDTreeCps *tree,
  double eps,
  // double (*randfun)(int),
  std::function<double (const int)> randfun,
  int *sample,
  int *sampleSize
) {
  IndexList *idx = new IndexList(N);
  int *neighbours = new int[N];
  double *weights = new double[N];
  double *dists = new double[N];

  idx->fill();

  while (idx->length() > 1) {
    int id1 = idx->draw();

    // We need to remove the unit first, so that it is not searching itself
    // in the tree search
    idx->erase(id1);
    tree->removeUnit(id1);

    // Find all neighbours
    int len = tree->findNeighbours(probabilities, weights, dists, neighbours, id1);

    bool included = randfun(id1) < probabilities[id1];
    double slag;

    if (included) {
      slag = probabilities[id1] - 1.0;
      probabilities[id1] = 1.0;

      *sampleSize += 1;
    } else {
      slag = probabilities[id1];
      // probabilities[id1] = 0.0;
    }

    double remweight = 1.0;

    // Loop through all found neighbours
    for (int i = 0; i < len && remweight > eps;) {
      // First we need to find how many neighbours exists on the same distance
      // Initialize totweight to the first neighbour, then search through
      // until the distance differs from this first neighbour
      double totweight = weights[neighbours[i]];

      int j = i + 1;
      for (; j < len; j++) {
        if (dists[neighbours[i]] < dists[neighbours[j]])
          break;

        totweight += weights[neighbours[j]];
      }

      // If we only found one nearest neighbour, we resolve this and continue
      if (j - i == 1) {
        int id2 = neighbours[i];

        // Do not use more than the remaining weight
        double temp = remweight >= totweight ? totweight : remweight;

        probabilities[id2] += temp * slag;
        decide(idx, tree, probabilities, sampleSize, id2, eps);

        i += 1;
        remweight -= temp;
        continue;
      }

      // If we found multiple nearest neighbours
      if (remweight >= totweight) {
        // The remaining weight is larger than the total weight of the nearest neighbours
        // Loop through the nearest neighbours and update their probabilities
        for (; i < j; i++) {
          int id2 = neighbours[i];
          probabilities[id2] += weights[id2] * slag;
          decide(idx, tree, probabilities, sampleSize, id2, eps);
        }

        remweight -= totweight;
      } else {
        // The remaining weight is smaller than the total weight of the nearest neighbours
        // We need to sort this list, smallest weights first
        std::sort(
          neighbours + i,
          neighbours + j,
          [weights](int a, int b) { return weights[a] < weights[b]; }
        );

        // Loop through all units, and update their weights
        // No unit can get more than a fair share
        for (; i < j; i++) {
          int id2 = neighbours[i];
          // Temp contains fair share
          double temp = remweight / (double)(j - i);
          // But we cannot update with more than the assigned weight
          if (weights[id2] < temp)
            temp = weights[id2];

          probabilities[id2] += temp * slag;
          decide(idx, tree, probabilities, sampleSize, id2, eps);
          remweight -= temp;
        }
      }
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (stduniform() < probabilities[id1]) {
      probabilities[id1] = 1.0;
      *sampleSize += 1;
    }
  }

  for (int i = 0, j = 0; i < N && j < *sampleSize; i++) {
    if (probabilities[i] >= 1.0 - eps) {
      sample[j] = i + 1;
      j += 1;
    }
  }

  delete[] neighbours;
  delete[] weights;
  delete[] dists;
  delete idx;

  return;
}

// [[Rcpp::export(.scps_cpp)]]
Rcpp::IntegerVector scps_cpp(
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

  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
  }

  std::function<double (const int)> randfun = [](int i) {
    return stduniform();
  };

  scps_internal(
    xx,
    N,
    probabilities,
    tree,
    eps,
    randfun,
    sample,
    &sampleSize
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete[] sample;
  delete tree;

  return svec;
}

// [[Rcpp::export(.scps_coord_cpp)]]
Rcpp::IntegerVector scps_coord_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  Rcpp::NumericVector &random,
  int bucketSize,
  int method,
  double eps
) {
  int N = x.ncol();
  double *xx = REAL(x);

  double *probabilities = new double[N];
  int *sample = new int[N];
  int sampleSize = 0;

  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
  }

  std::function<double (const int)> randfun = [&random](int i) {
    return random[i];
  };

  scps_internal(
    xx,
    N,
    probabilities,
    tree,
    eps,
    randfun,
    sample,
    &sampleSize
  );

  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] probabilities;
  delete[] sample;
  delete tree;

  return svec;
}
