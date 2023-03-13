#include <Rcpp.h>
#include "uniform.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))

// [[Rcpp::export(.spm_cpp)]]
Rcpp::IntegerVector spm_cpp(
  Rcpp::NumericVector &prob,
  double eps
) {
  int N = prob.length();
  int sampleSize = 0;
  int *sample = new int[N];
  int *idx = new int[N];
  double *probabilities = new double[N];

  for (int i = 0; i < N; i++) {
    idx[i] = i;
    probabilities[i] = prob[i];
  }

  int i = 0;
  while (i < N - 1) {
    int id1 = idx[i];
    int id2 = idx[i + 1];

    double p1 = probabilities[id1];
    double p2 = probabilities[id2];
    double psum = p1 + p2;

    if (psum > 1.0) {
      if (1.0 - p2 > stduniform() * (2.0 - psum)) {
        probabilities[id1] = 1.0;
        probabilities[id2] = psum - 1.0;
      } else {
        probabilities[id1] = psum - 1.0;
        probabilities[id2] = 1.0;
      }
    } else {
      if (p2 > stduniform() * psum) {
        probabilities[id1] = 0.0;
        probabilities[id2] = psum;
      } else {
        probabilities[id1] = psum;
        probabilities[id2] = 0.0;
      }
    }

    int decide = 0;

    if (pclose(probabilities[id1], eps)) {
      decide += 1;

      if (pbig(probabilities[id1], eps)) {
        sample[sampleSize] = id1 + 1;
        sampleSize += 1;
      }
    }

    if (pclose(probabilities[id2], eps)) {
      decide += 2;

      if (pbig(probabilities[id2], eps)) {
        sample[sampleSize] = id2 + 1;
        sampleSize += 1;
      }
    }

    switch(decide) {
    case 3: // Both were decide, we can just move up 2
      i += 2;
      break;
    case 2: // Only second was decided, we switch and move up 1
      idx[i + 1] = idx[i];
      i += 1;
      break;
    case 1: // Only first was decided, we move up 1
      i += 1;
      break;
    }
  }

  if (i == N - 1) {
    int id1 = idx[i];
    if (stduniform() < probabilities[id1]) {
      sample[sampleSize] = id1 + 1;
      sampleSize += 1;
    }
  }

  std::sort(sample, sample + sampleSize);
  Rcpp::IntegerVector svec(sample, sample + sampleSize);

  delete[] sample;
  delete[] idx;
  delete[] probabilities;

  return svec;
}
