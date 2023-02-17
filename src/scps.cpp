#include <algorithm>
#include <unordered_set>
#include <Rcpp.h>
#include "kdtree-cps.h"
#include "uniform.h"
#include "unorderedset-draw.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define intuniform(N) ((int)((double)N * stduniform()))
// #define pclose(p, eps) (p <= eps || p >= 1.0 - eps)

void decide(std::unordered_set<int> &idx, KDTreeCps *tree, double *probs, int *sample, int *nsample, int *unresolved, int uid, double eps) {
  if (probs[uid] <= eps) {
    idx.erase(uid);
    tree->removeUnit(uid);
    *unresolved -= 1;
  } else if (probs[uid] >= 1.0 - eps) {
    idx.erase(uid);
    tree->removeUnit(uid);
    *unresolved -= 1;
    sample[*nsample] = uid + 1;
    *nsample += 1;
  }
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
  int unresolvedObjects = N;
  double *xx = REAL(x);

  std::unordered_set<int> idx(N);
  int *neighbours = new int[N];
  double *probabilities = new double[N];
  double *weights = new double[N];
  double *dists = new double[N];
  int *sample = new int[N];
  int nsample = 0;

  KDTreeCps *tree = new KDTreeCps(xx, N, x.nrow(), bucketSize, method);
  tree->init();

  for (int i = 0; i < N; i++) {
    probabilities[i] = prob[i];
    idx.insert(i);
  }

  while (unresolvedObjects > 0) {
    int idx1 = unorderedsetDraw(idx, N);

    double slag = probabilities[idx1];
    int included = 0;

    if (stduniform() < probabilities[idx1]) {
      included = 1;
      sample[nsample] = idx1 + 1;
      nsample += 1;
      slag -= 1.0;
    }

    idx.erase(idx1);
    tree->removeUnit(idx1);
    unresolvedObjects -= 1;

    if (unresolvedObjects == 0)
      break;

    int len = tree->findNeighbours(probabilities, weights, dists, neighbours, idx1);
    double weight = 1.0;

    for (int i = 0; i < len && weight > eps;) {
      int j = i + 1;
      double totweight = weights[neighbours[i]];
      for (; j < len; j++) {
        if (dists[neighbours[i]] < dists[neighbours[j]])
          break;

        totweight += weights[neighbours[j]];
      }

      if (j - i == 1) {
        int idx2 = neighbours[i];
        double temp = weight >= totweight ? totweight : weight;
        probabilities[idx2] += temp * slag;
        decide(idx, tree, probabilities, sample, &nsample, &unresolvedObjects, idx2, eps);
        i += 1;
        weight -= temp;
        continue;
      }

      if (weight >= totweight) {
        for (; i < j; i++) {
          int idx2 = neighbours[i];
          probabilities[idx2] += weights[idx2] * slag;
          decide(idx, tree, probabilities, sample, &nsample, &unresolvedObjects, idx2, eps);
        }

        weight -= totweight;
      } else {
        std::sort(
          neighbours + i,
          neighbours + j,
          [weights](int a, int b) { return weights[a] < weights[b]; }
        );

        for (; i < j; i++) {
          int idx2 = neighbours[i];
          double temp = weight / (double)(j - i);
          if (weights[idx2] < temp) {
            temp = weights[idx2];
          }

          probabilities[idx2] += temp * slag;
          decide(idx, tree, probabilities, sample, &nsample, &unresolvedObjects, idx2, eps);
          weight -= temp;
        }
      }
    }
  }

  std::sort(sample, sample + nsample);
  Rcpp::IntegerVector sret(sample, sample + nsample);
  delete[] neighbours;
  delete[] probabilities;
  delete[] weights;
  delete[] dists;
  delete tree;

  return sret;
}
