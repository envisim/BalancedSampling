#include "cps-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define psmall(p, eps) ((p) <= (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))

CpsMethod intToCpsMethod(const int i) {
  if (1 <= i && i <= 3)
    return static_cast<CpsMethod>(i);

  std::invalid_argument("cps-method does not exist");
  return CpsMethod::err;
}

// DIRECT
Cps::Cps(
    const double* t_probabilities,
    double* xx,
    double* t_random,
    const int t_N,
    const int p,
    const CpsMethod cpsMethod,
    const int bucketSize,
    const int method,
    const double t_eps
  ) {
  set_indirect = false;

  N = t_N;
  eps = t_eps;

  sample = new int[N];
  probabilities = new double[N];
  randomValues = t_random;

  tree = new KDTreeCps(xx, N, p, bucketSize, method);
  tree->init();

  idx = new IndexList(N);

  if (N > 0) {
    // Decrement done before evaluating the loop, so it begins on N - 1
    for (int i = N; i-- > 0; ) {
      probabilities[i] = t_probabilities[i];
      idx->set(i);

      if (pclose(probabilities[i], eps)) {
        eraseUnit(i);

        if (pbig(probabilities[i], eps))
          addUnitToSample(i);
      }
    }
  }

  init(cpsMethod);
}

// INDIRECT
Cps::Cps(
    double* t_probabilities,
    KDTreeCps* t_tree,
    IndexList* t_idx,
    double* t_random,
    const int t_N,
    const CpsMethod cpsMethod,
    const double t_eps
  ) {
  set_indirect = true;

  probabilities = t_probabilities;
  tree = t_tree;
  idx = t_idx;
  randomValues = t_random;
  N = t_N;
  eps = t_eps;

  sample = new int[N];
  init(cpsMethod);
}

void Cps::init(const CpsMethod cpsMethod) {
  neighbours = new int[N];
  distances = new double[N];
  weights = new double[N];

  if (randomValues == nullptr) {
    _random = &Cps::random_std;
  } else {
    _random = &Cps::random_arr;
  }

  switch(cpsMethod) {
  case CpsMethod::lcps:
    neighbours2 = new int[N];
    _draw = &Cps::draw_lcps;
    break;

  case CpsMethod::scps:
    _draw = &Cps::draw_scps;
    break;

  case CpsMethod::scpscoord:
    _draw = &Cps::draw_scpscoord;
    if (randomValues == nullptr)
      std::invalid_argument("random not set for scpscoord");
    break;

  default:
    std::invalid_argument("cpsMethod does not exist");
    break;
  }

  set_draw = true;

  return;
}

void Cps::addUnitToSample(const int id) {
  sample[sampleSize] = id + 1;
  sampleSize += 1;
  return;
}

void Cps::eraseUnit(const int id) {
  idx->erase(id);
  tree->removeUnit(id);
  return;
}

void Cps::decideUnit(const int id) {
  if (pclose(probabilities[id], eps)) {
    eraseUnit(id);

    if (pbig(probabilities[id], eps))
      addUnitToSample(id);
  }

  return;
}

int Cps::draw_lcps() {
  // Take care of edge cases
  if (idx->length() <= 1) {
    if (idx->length() == 1)
      return idx->get(0);
    if (idx->length() < 1)
      std::range_error("trying to find index in empty list");
  }

  int arrlen = 0;
  double mindist = DBL_MAX;

  // Loop through all remaining units.
  // Put the smallest distances in idxarr
  for (int i = 0; i < idx->length(); i++) {
    int len = tree->findNeighbours(probabilities, weights, distances, neighbours, idx->get(i));
    double dist = distances[neighbours[len - 1]];

    if (dist < mindist) {
      neighbours2[0] = i;
      mindist = dist;
      arrlen = 1;
    } else if (dist == mindist) {
      neighbours2[arrlen] = i;
      arrlen += 1;
    }
  }

  // Choose randomly from the units in idxarr
  int idxi = intuniform(arrlen);
  return idx->get(idxi);
}

int Cps::draw_scps() {
  return idx->draw();
}

int Cps::draw_scpscoord() {
  while(!idx->exists(histn))
    histn += 1;

  int tunit = histn;
  histn += 1;
  return tunit;
}

double Cps::random_std(const int id) {
  return stduniform();
}

double Cps::random_arr(const int id) {
  return randomValues[id];
}

void Cps::_run() {
  while (idx->length() > 1) {
    int id1 = draw();

    // We need to remove the unit first, so that it is not searching itself
    // in the tree search
    eraseUnit(id1);

    // Find all neighbours
    int len = tree->findNeighbours(probabilities, weights, distances, neighbours, id1);

    double slag = probabilities[id1];

    if (random(id1) < probabilities[id1]) {
      slag -= 1.0;
      addUnitToSample(id1);
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
        if (distances[neighbours[i]] < distances[neighbours[j]])
          break;

        totweight += weights[neighbours[j]];
      }

      // If we only found one nearest neighbour, we resolve this and continue
      if (j - i == 1) {
        int id2 = neighbours[i];

        // Do not use more than the remaining weight
        double temp = remweight >= totweight ? totweight : remweight;

        probabilities[id2] += temp * slag;
        decideUnit(id2);

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
          decideUnit(id2);
        }

        remweight -= totweight;
      } else {
        // The remaining weight is smaller than the total weight of the nearest neighbours
        // We need to sort this list, smallest weights first
        std::sort(
          neighbours + i,
          neighbours + j,
          [this](int a, int b) { return weights[a] < weights[b]; }
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
          decideUnit(id2);
          remweight -= temp;
        }
      }
    }
  }

  // Sort out any remaining lone unit
  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (random(id1) < probabilities[id1])
      addUnitToSample(id1);
  }

  std::sort(sample, sample + sampleSize);
  return;
}
