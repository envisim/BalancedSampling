#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))

void lpm_internal(
  KDTree* tree,
  IndexList* idx,
  double* probabilities,
  const int N,
  const double eps,
  int* sample,
  int* sampleSize,
  std::function<void (int*)> unitfun
) {
  int *pair = new int[2];

  while (idx->length() > 1) {
    unitfun(pair);
    int id1 = pair[0];
    int id2 = pair[1];

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

    if (pclose(probabilities[id1], eps)) {
      idx->erase(id1);
      tree->removeUnit(id1);

      if (pbig(probabilities[id1], eps)) {
        sample[*sampleSize] = id1 + 1;
        *sampleSize += 1;
      }
    }

    if (pclose(probabilities[id2], eps)) {
      idx->erase(id2);
      tree->removeUnit(id2);

      if (pbig(probabilities[id2], eps)) {
        sample[*sampleSize] = id2 + 1;
        *sampleSize += 1;
      }
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (stduniform() < probabilities[id1]) {
      sample[*sampleSize] = id1 + 1;
      *sampleSize += 1;
    }
  }

  delete[] pair;

  std::sort(sample, sample + *sampleSize);
  return;
}

void lpm_int_internal(
  KDTree* tree,
  IndexList* idx,
  int* probabilities,
  const int N,
  int* sample,
  int* sampleSize,
  std::function<void (int*)> unitfun
) {
  int *pair = new int[2];

  while (idx->length() > 1) {
    unitfun(pair);
    int id1 = pair[0];
    int id2 = pair[1];

    int p1 = probabilities[id1];
    int p2 = probabilities[id2];
    int psum = p1 + p2;

    if (psum > N) {
      if (N - p2 > intuniform(2 * N - psum)) {
        probabilities[id1] = N;
        probabilities[id2] = psum - N;
      } else {
        probabilities[id1] = psum - N;
        probabilities[id2] = N;
      }
    } else {
      if (p2 > intuniform(psum)) {
        probabilities[id1] = 0;
        probabilities[id2] = psum;
      } else {
        probabilities[id1] = psum;
        probabilities[id2] = 0;
      }
    }

    if (probabilities[id1] == 0 || probabilities[id1] == N) {
      idx->erase(id1);
      tree->removeUnit(id1);

      if (probabilities[id1] == N) {
        sample[*sampleSize] = id1 + 1;
        *sampleSize += 1;
      }
    }

    if (probabilities[id2] == 0 || probabilities[id2] == N) {
      idx->erase(id2);
      tree->removeUnit(id2);

      if (probabilities[id2] == N) {
        sample[*sampleSize] = id2 + 1;
        *sampleSize += 1;
      }
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (intuniform(N) < probabilities[id1]) {
      sample[*sampleSize] = id1 + 1;
      *sampleSize += 1;
    }
  }

  delete[] pair;

  std::sort(sample, sample + *sampleSize);
  return;
}
