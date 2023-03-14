#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))

// DOUBLE
Lpm::Lpm(
  const double* t_probabilities,
  double* xx,
  const int t_N,
  const int p,
  const int lpMethod,
  const int bucketSize,
  const int method,
  const double t_eps
) {
  set_indirect = false;

  N = t_N;
  eps = t_eps;

  probabilities = new double[N];

  tree = new KDTree(xx, N, p, bucketSize, method);
  tree->init();

  idx = new IndexList(N);

  for (int i = 0; i < N; i++) {
    probabilities[i] = t_probabilities[i];
    idx->set(i);
  }

  init(lpMethod, false);
}
// INT
Lpm::Lpm(
  const int t_pn,
  double* xx,
  const int t_N,
  const int p,
  const int lpMethod,
  const int bucketSize,
  const int method
) {
  set_indirect = false;

  N = t_N;

  iprobabilities = new int[N];

  tree = new KDTree(xx, N, p, bucketSize, method);
  tree->init();

  idx = new IndexList(N);

  for (int i = 0; i < N; i++) {
    iprobabilities[i] = t_pn;
    idx->set(i);
  }

  init(lpMethod, true);
}
// DOUBLE SET
Lpm::Lpm(
  double* t_probabilities,
  KDTree* t_tree,
  IndexList* t_idx,
  const int t_N,
  const int lpMethod,
  const double t_eps
) {
  set_indirect = true;

  probabilities = t_probabilities;
  tree = t_tree;
  idx = t_idx;
  N = t_N;
  eps = t_eps;

  init(lpMethod, false);
}
// INT SET
Lpm::Lpm(
  int* t_probabilities,
  KDTree* t_tree,
  IndexList* t_idx,
  const int t_N,
  const int lpMethod
) {
  set_indirect = true;

  iprobabilities = t_probabilities;
  tree = t_tree;
  idx = t_idx;
  N = t_N;

  init(lpMethod, true);
}

void Lpm::init(const int lpMethod, const bool isInt) {
  sample = new int[N];
  neighbours = new int[N];

  // sampleSize = 0;

  if (isInt) {
    _run = &Lpm::run_int;
  } else {
    _run = &Lpm::run_double;
  }

  set_run = true;

  if (lpMethod == 1) {
    neighbours2 = new int[N];
    _draw = &Lpm::draw_lpm1;
  } else if (lpMethod == 3) {
    neighbours2 = new int[N];
    history = new int[N];
    _draw = &Lpm::draw_lpm1search;
  } else {
    // lpMethod == 2
    _draw = &Lpm::draw_lpm2;
  }

  set_draw = true;

  return;
}

void Lpm::addUnitToSample(const int id) {
  sample[sampleSize] = id + 1;
  sampleSize += 1;
  return;
}

void Lpm::eraseUnit(const int id) {
  idx->erase(id);
  tree->removeUnit(id);
  return;
}

void Lpm::draw_lpm1(int* pair) {
  while (true) {
    pair[0] = idx->draw();
    int len = tree->findNeighbour(neighbours, N, pair[0]);

    for (int i = 0; i < len;) {
      int tlen = tree->findNeighbour(neighbours2, N, neighbours[i]);
      bool found = false;

      for (int j = 0; j < tlen; j++) {
        if (neighbours2[j] == pair[0]) {
          found = true;
          break;
        }
      }

      if (found) {
        i += 1;
      } else {
        len -= 1;
        neighbours[i] = neighbours[len];
      }
    }

    if (len > 0) {
      pair[1] = neighbours[intuniform(len)];
      return;
    }
  }
}

void Lpm::draw_lpm2(int* pair) {
  pair[0] = idx->draw();
  int len = tree->findNeighbour(neighbours, N, pair[0]);
  pair[1] = neighbours[intuniform(len)];
  return;
}

void Lpm::draw_lpm1search(int* pair) {
  // Go back in the history and remove units that does not exist
  while (histn > 0) {
    if (idx->exists(history[histn - 1])) {
      break;
    }

    histn -= 1;
  }

  // If there is no history, we draw a unit at random
  if (histn == 0) {
    history[0] = idx->draw();
    histn = 1;
  }

  while (true) {
    // Set the first unit to the last in history
    pair[0] = history[histn - 1];
    // Find this units nearest neighbours
    int len = tree->findNeighbour(neighbours, N, pair[0]);
    int len_copy = len;

    // Go through all nearest neighbours
    for (int i = 0; i < len;) {
      // Find the neighbours nearest neighbours
      int tlen = tree->findNeighbour(neighbours2, N, neighbours[i]);
      bool found = false;

      // Check if any of these are the history-unit
      for (int j = 0; j < tlen; j++) {
        if (neighbours2[j] == pair[0]) {
          found = true;
          break;
        }
      }

      // If the history-unit exists among the nearest neighbours, we continue
      // to see if any other of the history-units neighbours also are mutual.
      // Otherwise, the history-unit is not among the nearest neighbours,
      // we swap places and continue the search.
      if (found) {
        i += 1;
      } else {
        len -= 1;
        if (i != len) {
          int temp = neighbours[i];
          neighbours[i] = neighbours[len];
          neighbours[len] = temp;
        }
      }
    }

    // If we found one or more mutual neighbours, we select one at random
    if (len > 0) {
      pair[1] = neighbours[intuniform(len)];
      return;
    }

    // If we come here, no mutual neighbours exist

    // We might need to clear the history if the search has been going on for
    // too long. This can probably? happen if there is a long history, and
    // updates has affected previous units.
    if (histn == N)
      histn = 0;

    // We select a unit at random to become the next history unit, and traverse
    // one step further.
    history[histn] = neighbours[intuniform(len_copy)];
    histn += 1;
  }
}

void Lpm::run_double() {
  int *pair = new int[2];

  while (idx->length() > 1) {
    draw(pair);
    int id1 = pair[0];
    int id2 = pair[1];

    double* p1 = probabilities + id1;
    double* p2 = probabilities + id2;
    double psum = *p1 + *p2;

    if (psum > 1.0) {
      if (1.0 - *p2 > stduniform() * (2.0 - psum)) {
        *p1 = 1.0;
        *p2 = psum - 1.0;
      } else {
        *p1 = psum - 1.0;
        *p2 = 1.0;
      }
    } else {
      if (*p2 > stduniform() * psum) {
        *p1 = 0.0;
        *p2 = psum;
      } else {
        *p1 = psum;
        *p2 = 0.0;
      }
    }

    if (pclose(*p1, eps)) {
      eraseUnit(id1);

      if (pbig(*p1, eps)) {
        addUnitToSample(id1);
      }
    }

    if (pclose(*p2, eps)) {
      eraseUnit(id2);

      if (pbig(*p2, eps)) {
        addUnitToSample(id2);
      }
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (stduniform() < probabilities[id1]) {
      addUnitToSample(id1);
    }
  }

  delete[] pair;

  std::sort(sample, sample + sampleSize);
  return;
}

void Lpm::run_int() {
  int *pair = new int[2];

  while (idx->length() > 1) {
    draw(pair);
    int id1 = pair[0];
    int id2 = pair[1];

    int* p1 = iprobabilities + id1;
    int* p2 = iprobabilities + id2;
    int psum = *p1 + *p2;

    if (psum > N) {
      if (N - *p2 > intuniform(2 * N - psum)) {
        *p1 = N;
        *p2 = psum - N;
      } else {
        *p1 = psum - N;
        *p2 = N;
      }
    } else {
      if (*p2 > intuniform(psum)) {
        *p1 = 0;
        *p2 = psum;
      } else {
        *p1 = psum;
        *p2 = 0;
      }
    }

    if (*p1 == 0 || *p1 == N) {
      eraseUnit(id1);

      if (*p1 == N) {
        addUnitToSample(id1);
      }
    }

    if (*p2 == 0 || *p2 == N) {
      eraseUnit(id2);

      if (*p2 == N) {
        addUnitToSample(id2);
      }
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);
    if (intuniform(N) < iprobabilities[id1]) {
      addUnitToSample(id1);
    }
  }

  delete[] pair;

  std::sort(sample, sample + sampleSize);
  return;
}

