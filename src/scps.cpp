#include <algorithm>
#include <Rcpp.h>
#include "kdtree-cps.h"
#include "uniform.h"
#include "index-list.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define intuniform(N) ((int)((double)N * stduniform()))
// #define pclose(p, eps) (p <= eps || p >= 1.0 - eps)

void decide(
  IndexList *idx,
  KDTreeCps *tree,
  double *probs,
  int *sample,
  int *nsample,
  int uid,
  double eps
) {
  if (probs[uid] <= eps) {
    idx->erase(uid);
    tree->removeUnit(uid);
  } else if (probs[uid] >= 1.0 - eps) {
    idx->erase(uid);
    tree->removeUnit(uid);
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
  double *xx = REAL(x);

  IndexList *idx = new IndexList(N);
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
    idx->set(i);
  }

  while (idx->length() > 0) {
    int id1 = idx->draw();

    double slag = probabilities[id1];
    int included = 0;

    if (stduniform() < probabilities[id1]) {
      included = 1;
      sample[nsample] = id1 + 1;
      nsample += 1;
      slag -= 1.0;
    }

    idx->erase(id1);
    tree->removeUnit(id1);

    if (idx->length() == 0)
      break;

    int len = tree->findNeighbours(probabilities, weights, dists, neighbours, id1);
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
        int id2 = neighbours[i];
        double temp = weight >= totweight ? totweight : weight;
        probabilities[id2] += temp * slag;
        decide(idx, tree, probabilities, sample, &nsample, id2, eps);
        i += 1;
        weight -= temp;
        continue;
      }

      if (weight >= totweight) {
        for (; i < j; i++) {
          int id2 = neighbours[i];
          probabilities[id2] += weights[id2] * slag;
          decide(idx, tree, probabilities, sample, &nsample, id2, eps);
        }

        weight -= totweight;
      } else {
        std::sort(
          neighbours + i,
          neighbours + j,
          [weights](int a, int b) { return weights[a] < weights[b]; }
        );

        for (; i < j; i++) {
          int id2 = neighbours[i];
          double temp = weight / (double)(j - i);
          if (weights[id2] < temp) {
            temp = weights[id2];
          }

          probabilities[id2] += temp * slag;
          decide(idx, tree, probabilities, sample, &nsample, id2, eps);
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
  delete idx;
  delete tree;

  return sret;
}
