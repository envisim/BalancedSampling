#include "kdtree-cps.h"

int KDTreeCps::findNeighbours(
  const double *probabilities,
  double *weights,
  double *dists,
  int *neighbours,
  // const int maxsize,
  const int idx
) {
  if (top == nullptr)
    return 0;

  // Ready the pointer to the deciding unit
  double *unit = data + idx * p;
  // Ready the size of the neighbours-vector
  int size = 0;
  // mindist[0] = DBL_MAX;
  // Ready the total weight sum
  double weightsum = 0.0;

  findNeighboursSearch(probabilities, weights, dists, neighbours, &size, top, &weightsum, idx, unit);

  return size;
}

void KDTreeCps::findNeighboursSearch(
  const double *probabilities,
  double *weights,
  double *dists,
  int *neighbours,
  // const int maxsize,
  int *size,
  KDNode *node,
  double *weightsum,
  const int idx,
  double *unit
) {
  // If the current node is a terminal node, we investigate within
  if (node->isTerminal()) {
    findNeighboursInNode(probabilities, weights, dists, neighbours, size, node, weightsum, idx, unit);
    return;
  }

  // Decide primary search direction
  double dist = *(unit + node->split) - node->value;
  KDNode *nextnode = dist <= 0.0 ? node->cleft : node->cright;

  // Search within this node
  findNeighboursSearch(probabilities, weights, dists, neighbours, size, nextnode, weightsum, idx, unit);

  // <= to be precise, but maybe we can ignore this...
  // if (square(data[idx * p + node->split] - node->value) < *mindist) {
  if (*weightsum < 1.0 || dist * dist < dists[neighbours[*size - 1]]) {
    findNeighboursSearch(probabilities, weights, dists, neighbours, size, nextnode->getSibling(), weightsum, idx, unit);
  }
}

void KDTreeCps::findNeighboursInNode(
  const double *probabilities,
  double *weights,
  double *dists,
  int *neighbours,
  // const int maxsize,
  int *size,
  KDNode *node,
  double *weightsum,
  const int idx,
  double *unit
) {
  int nunits = node->getSize();

  // If there are no units within this node, we need not to bother
  if (nunits == 0)
    return;

  int tempsize = *size;
  double distmin = DBL_MAX;

  // Search through all units in the node, and decide the distances and weights
  for (int i = 0; i < nunits; i++) {
    int idx2 = node->units[i];
    // Can't be same unit, always decided
    // if (idx2 == idx) // Same unit
    //   continue;

    neighbours[tempsize] = idx2;
    tempsize += 1;

    dists[idx2] = distance(unit, data + idx2 * p);
    weights[idx2] = probabilities[idx] + probabilities[idx2] <= 1.0
      ? probabilities[idx2] / (1.0 - probabilities[idx])
      : (1.0 - probabilities[idx2]) / probabilities[idx];

    if (dists[idx2] < distmin)
      distmin = dists[idx2];
  }

  // If weightsum is full and the current border dist is less than the node's
  // distmin, then we have nothing here to do
  if (*weightsum >= 1.0 && *size > 0 && dists[neighbours[*size - 1]] < distmin)
    return;

  // Otherwise, we have four possibilities
  //   (A)
  // - Either there is no units in neighbours, then we need to sort and add all
  //     the newly acquired units
  //   (B)
  // - Or, the smallest distance in the node is less than the current smallest,
  //     and we also need to sort and add all the newly acquired units, together
  //     with the previously aquired
  //   (C)
  // - Or, the smallest newly acquired unit is larger than the current largest,
  //     while we already know that the weightsum is < 1.0, so we need to only
  //     sort the newly acquired units
  //   (D)
  // - Finally, the smallest newly acquired unit is lying somewhere within the
  //     current units, and we need to find a place to start the sorting on.

  // in neighbours, we have two groups:
  // - sorted current units [0, size)
  // - unsorted newly acquired units [size, tempsize)
  // i decides where to start sorting
  int i;
  if (*size == 0 || distmin < dists[neighbours[0]]) {
    // (A) and (B), reset weightsum, we need to sort [0, tempsize)
    i = 0;
    *weightsum = 0.0;
  } else if (dists[neighbours[*size - 1]] <= distmin) {
    // (C), keep all old, sort [size, tempsize)
    i = *size;
    // *weightsum = *weightsum;
  } else {
    // (D), go backwards in [0, size), and remove weights for all units that
    //      need sorting
    for (i = *size - 1; i >= 0; i--) {
      if (dists[neighbours[i]] < distmin)
        break;

      *weightsum -= weights[neighbours[i]];
    }

    i += 1;
  }

  // Sort the range [i, tempsize)
  std::sort(
    neighbours + i,
    neighbours + tempsize,
    [dists](int a, int b) { return dists[a] < dists[b]; }
  );

  // Add unit i to the list
  *weightsum += weights[neighbours[i]];
  // i += 1;

  // Add any following unit, stopping when
  // - weightsum >= 1.0, and
  // - the previous distance is smaller than the current distance
  for (i = i + 1; i < tempsize; i++) {
    if (*weightsum >= 1.0 && dists[neighbours[i - 1]] < dists[neighbours[i]])
      break;

    *weightsum += weights[neighbours[i]];
  }

  *size = i;
  return;
}
