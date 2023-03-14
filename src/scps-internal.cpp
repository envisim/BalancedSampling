#include "scps-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

void scpsDecide(
  ScpsDecideProps* props,
  const int uid
) {
  if (props->probs[uid] <= props->eps) {
    props->idx->erase(uid);
    props->tree->removeUnit(uid);
  } else if (props->probs[uid] >= 1.0 - props->eps) {
    props->idx->erase(uid);
    props->tree->removeUnit(uid);
    props->sample[*(props->sampleSize)] = uid + 1;
    *(props->sampleSize) += 1;
  }
}

void scps_internal(
  KDTreeCps* tree,
  IndexList* idx,
  double* probabilities,
  const int N,
  const double eps,
  int* sample,
  int* sampleSize,
  std::function<double (const int)> randfun,
  std::function<int ()> unitfun
) {
  // Initialize arrays needed for neighbour search
  int *neighbours = new int[N];
  double *weights = new double[N];
  double *dists = new double[N];

  // Initialize decideProps, for use in scpsDecideUnits
  ScpsDecideProps decideProps(tree, idx, probabilities, eps, sample, sampleSize);

  // Loop through this until at most one unit exists
  while (idx->length() > 1) {
    int id1 = unitfun();

    // We need to remove the unit first, so that it is not searching itself
    // in the tree search
    idx->erase(id1);
    tree->removeUnit(id1);

    // Find all neighbours
    int len = tree->findNeighbours(probabilities, weights, dists, neighbours, id1);

    double slag = probabilities[id1];

    if (randfun(id1) < probabilities[id1]) {
      slag -= 1.0;
      sample[*sampleSize] = id1 + 1;
      *sampleSize += 1;
    }

    // The weight that remains to be put out to the neighbours
    double remweight = 1.0;

    // Loop through all found neighbours
    // The loop is conducted so that we take equal distance neighbours together
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

  // Sort out any remaining lone unit
  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (randfun(id1) < probabilities[id1]) {
      sample[*sampleSize] = id1 + 1;
      *sampleSize += 1;
    }
  }

  delete[] neighbours;
  delete[] weights;
  delete[] dists;

  // Sort the sample list
  std::sort(sample, sample + *sampleSize);
  return;
}
