#include "lpm2-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define intuniform(N) ((int)((double)N * stduniform()))
#define pbig(p, eps) (p >= 1.0 - eps)
#define pclose(p, eps) (p <= eps || p >= 1.0 - eps)

void lpm2_internal(
  KDTree *tree,
  IndexList *idx,
  double *probabilities,
  int *neighbours,
  const int N,
  const double eps,
  int *sample,
  int *sampleSize
) {
  while (idx->length() > 1) {
    int id1 = idx->draw();

    int len = tree->findNeighbour(neighbours, N, id1);
    int id2 = len == 1 ? neighbours[0] : neighbours[intuniform(len)];

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

  std::sort(sample, sample + *sampleSize);
  return;
}
