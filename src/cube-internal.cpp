#include "cube-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))
#define mati(r, c, p) ((r) * (p) + (c))

void cubeUpdatePhase(
  double *bmat,
  double *probabilities,
  double *uvec,
  int *index,
  IndexList* idx,
  int *sample,
  int *sampleSize,
  const double eps,
  const int maxSize
) {
  rref(bmat, maxSize - 1, maxSize);

  double lambda1 = DBL_MAX;
  double lambda2 = DBL_MAX;

  for (int i = 0; i < maxSize; i++) {
    if (i == maxSize - 1) {
      uvec[i] = 1.0;
    } else {
      uvec[i] = -bmat[mati(i, maxSize - 1, maxSize)];
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

  for (int i = 0; i < maxSize; i++) {
    int id = index[i];
    probabilities[id] += lambda * uvec[i];

    if (pclose(probabilities[id], eps)) {
      idx->erase(id);

      if (pbig(probabilities[id], eps)) {
        sample[*sampleSize] = id + 1;
        *sampleSize += 1;
      }
    }
  }

  return;
}

void cubeFlightPhase(
  const double *amat, // To copy from
  double *probabilities, // To set
  IndexList* idx,
  KDTree* tree,
  int *sample,
  int *sampleSize,
  const double eps,
  const int p,
  const int N,
  int *index,
  double *distances,
  double *uvec,
  double *bmat
) {
  // Cases:
  // - choose from tree and all
  // - choose from list and all

  int maxSize = p + 1;

  while (idx->length() >= maxSize) {
    if (tree == nullptr) {
      for (int i = 0; i < maxSize; i++)
        index[i] = idx->get(i);
    } else {
      int id = idx->draw();
      index[0] = id;
      int len = tree->findNeighboursN(distances, index + 1, maxSize - 1, id);
      len += 1;

      if (len > maxSize) {
        // We have equals at the end
        int i = len - 1;
        while (i > 1) {
          if (distances[index[i]] > distances[index[i - 1]])
            break;

          i -= 1;
        }

        for (; i < maxSize; i++) {
          int j = intuniform(len - i) + i;
          if (i == j)
            continue;
          int temp = index[i];
          index[i] == index[j];
          index[j] = temp;
        }
      }
    }

    for (int i = 0; i < maxSize; i++) {
      for (int k = 0; k < maxSize - 1; k++) {
        bmat[mati(k, i, maxSize)] = amat[mati(k, index[i], N)];
      }
    }

    cubeUpdatePhase(
      bmat,probabilities,
      uvec,
      index,
      idx,
      sample,
      sampleSize,
      eps,
      maxSize
    );

    if (tree != nullptr) {
      for (int i = 0; i < maxSize; i++) {
        if (idx->exists(index[i]))
          continue;

        tree->removeUnit(index[i]);
      }
    }
  }

  return;
}

void cubeLandingPhase(
  const double *amat,
  double *probabilities,
  IndexList* idx,
  // KDTree* tree,
  int *sample,
  int *sampleSize,
  const double eps,
  const int p,
  const int N,
  int *index,
  double *uvec,
  double *bmat
) {
  if (idx->length() >= p + 1)
    std::range_error("landingphase committed early");

  while (idx->length() > 1) {
    int maxSize = idx->length();

    for (int i = 0; i < maxSize; i++) {
      index[i] = idx->get(i);

      for (int k = 0; k < maxSize - 1; k++) {
        bmat[mati(k, i, maxSize)] = amat[mati(k, index[i], N)];
      }
    }

    cubeUpdatePhase(
      bmat,
      probabilities,
      uvec,
      index,
      idx,
      sample,
      sampleSize,
      eps,
      maxSize
    );

    // if (tree != nullptr) {
    //   for (int i = 0; i < maxSize; i++) {
    //     if (idx->exists(index[i]))
    //       continue;

    //     tree->removeUnit(index[i]);
    //   }
    // }
  }

  if (idx->length() == 1) {
    int id = idx->get(0);
    if (stduniform() < probabilities[id]) {
      sample[*sampleSize] = id + 1;
      *sampleSize += 1;
    }
  }

  return;
}

void cubeRunPhases(
  const double *amat,
  double *probabilities,
  IndexList* idx,
  KDTree* tree,
  int *sample,
  int *sampleSize,
  const double eps,
  const int p,
  const int N
) {
  int maxSize = p + 1;

  int *index = new int[N];
  double *distances = tree == nullptr ? nullptr : new double[N];
  double *uvec = new double[maxSize];
  double *bmat = new double[maxSize * (maxSize - 1)];

  cubeFlightPhase(
    amat,
    probabilities,
    idx,
    tree,
    sample,
    sampleSize,
    eps,
    p,
    N,
    index,
    distances,
    uvec,
    bmat
  );

  cubeLandingPhase(
    amat,
    probabilities,
    idx,
    // tree,
    sample,
    sampleSize,
    eps,
    p,
    N,
    index,
    uvec,
    bmat
  );

  delete[] index;
  delete[] distances;
  delete[] uvec;
  delete[] bmat;

  return;
}

