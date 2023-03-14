#include "scps-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

void scpsDecide(
  ScpsDecideProps *props,
  const int uid
) {
  if (props->probs[uid] <= props->eps) {
    props->idx->erase(uid);
    props->tree->removeUnit(uid);
  } else if (props->probs[uid] >= 1.0 - props->eps) {
    props->idx->erase(uid);
    props->tree->removeUnit(uid);
    *(props->sampleSize) += 1;
  }
}

void scps_internal(
  const double *xx,
  const int N,
  double *probabilities,
  KDTreeCps *tree,
  const double eps,
  std::function<double (const int)> randfun,
  std::function<int ()> unitfun,
  int *sample,
  int *sampleSize,
  IndexList *idx
) {
  int *neighbours = new int[N];
  double *weights = new double[N];
  double *dists = new double[N];

  ScpsDecideProps decideProps(idx, tree, probabilities, sampleSize, eps);

  while (idx->length() > 1) {
    // int id1 = idx->draw();
    int id1 = unitfun();

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
        scpsDecide(&decideProps, id2);

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
          scpsDecide(&decideProps, id2);
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
          scpsDecide(&decideProps, id2);
          remweight -= temp;
        }
      }
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (randfun(id1) < probabilities[id1]) {
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

  return;
}
