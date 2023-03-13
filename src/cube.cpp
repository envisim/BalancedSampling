#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include "cube-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))
#define mati(r, c, p) ((r) * (p) + (c))

// [[Rcpp::export(.cube_cpp)]]
Rcpp::IntegerVector cube_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  double eps
) {
  int N = x.nrow();
  int p0 = x.ncol();
  int p = p0;
  double *amat = new double[N * p];
  double *bmat = new double[N * p];
  double *probabilities = new double[N];
  int *index = new int[N];
  int *sample = new int[N];
  int sampleSize = 0;
  int decidedSize = 0;

  double *uvec = new double[N];

  for (int i = N - 1; i >= 0; i--) {
    if (pclose(prob[i], eps)) {
      if (pbig(prob[i], eps)) {
        sample[sampleSize] = i;
        sampleSize += 1;
      }

      decidedSize += 1;
      if (i != N - decidedSize) {
        index[i] = index[N - decidedSize];
      }

      continue;
    }

    index[i] = i;
    probabilities[i] = prob[i];

    for (int k = 0; k < p0; k++) {
      amat[mati(k, i, N)] = x(i, k) / prob[i];
    }
  }

  while (decidedSize < N - 1) {
    if (decidedSize + 1 > N - p) {
      p -= 1;
    }

    int remaining = N - decidedSize;

    for (int i = 0; i < remaining; i++) {
      for (int k = 0; k < p; k++) {
        bmat[mati(k, i, remaining)] = amat[mati(k, index[i], N)];
      }
    }

    rref(bmat, p, remaining);

    double lambda1 = DBL_MAX;
    double lambda2 = DBL_MAX;

    for (int i = 0; i < remaining; i++) {
      if (i < p) {
        uvec[i] = 0.0;
        for (int k = p; k < remaining; k++) {
          // uvec[i] -= bmat[k * p + i];
          uvec[i] -= bmat[mati(i, k, remaining)];
        }
      } else {
        uvec[i] = 1.0;
      }

      double lval1 = std::abs(probabilities[index[i]] / uvec[i]);
      double lval2 = std::abs((1.0 - probabilities[index[i]]) / uvec[i]);

      if (uvec[i] >= 0.0) {
        if (lambda1 > lval2)
          lambda1 = lval2;
        if (lambda2 > lval1)
          lambda2 = lval1;
      } else {
        if (lambda1 > lval1)
          lambda1 = lval1;
        if (lambda2 > lval2)
          lambda2 = lval2;
      }
    }

    double lambda = stduniform() * (lambda1 + lambda2) < lambda2 ?
      lambda1 : -lambda2;

    for (int i = 0; i < remaining; i++) {
      int id = index[i];
      probabilities[id] += lambda * uvec[i];

      if (pclose(probabilities[id], eps)) {
        decidedSize += 1;

        if (id != N - decidedSize) {
          index[i] = index[N - decidedSize];
        }

        if (pbig(probabilities[id], eps)) {
          sample[sampleSize] = id + 1;
          sampleSize += 1;
        }
      }
    }
  }

  if (decidedSize == N - 1) {
    int id = index[0];
    if (stduniform() < probabilities[id]) {
      sample[sampleSize] = id + 1;
      sampleSize += 1;
    }
  }

  std::sort(sample, sample + sampleSize);
  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] amat;
  delete[] bmat;
  delete[] probabilities;
  delete[] index;
  delete[] sample;
  delete[] uvec;

  return svec;
}

// [[Rcpp::export(.cube_fast_cpp)]]
Rcpp::IntegerVector cube_fast_cpp(
  Rcpp::NumericVector &prob,
  Rcpp::NumericMatrix &x,
  double eps
) {
  int N = x.nrow();
  int p = x.ncol();
  int sampleSize = 0;

  double *amat = new double[N * p];
  double *probabilities = new double[N];
  IndexList *idx = new IndexList(N);
  int *sample = new int[N];

  for (int i = 0; i < N; i++) {
    idx->set(i);
    probabilities[i] = prob[i];
  }

  for (int i = 0; i < N; i++) {
    if (pclose(prob[i], eps)) {
      if (pbig(prob[i], eps)) {
        sample[sampleSize] = i;
        sampleSize += 1;
      }

      idx->erase(i);

      continue;
    }

    for (int k = 0; k < p; k++) {
      amat[mati(k, i, N)] = x(i, k) / prob[i];
    }
  }

  idx->shuffle();

  cubeRunPhases(
    amat,
    probabilities,
    idx,
    nullptr,
    sample,
    &sampleSize,
    eps,
    p,
    N
  );

  std::sort(sample, sample + sampleSize);
  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] amat;
  delete[] probabilities;
  delete idx;
  delete[] sample;

  return svec;
}


