#include "lpm-internal.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

#define pbig(p, eps) ((p) >= 1.0 - (eps))
#define pclose(p, eps) ((p) <= (eps) || (p) >= 1.0 - (eps))

LpmMethod intToLpmMethod(const int i) {
  if (1 <= i && i <= 5)
    return static_cast<LpmMethod>(i);

  std::invalid_argument("lpm-method does not exist");
  return LpmMethod::err;
}

// DOUBLE DIRECT
Lpm::Lpm(
  const double* t_probabilities,
  double* xx,
  const int t_N,
  const int p,
  const LpmMethod lpMethod,
  const int bucketSize,
  const int method,
  const double t_eps
) {
  set_indirect = false;

  N = t_N;
  eps = t_eps;

  sample = new int[N];
  probabilities = new double[N];

  tree = new KDTree(xx, N, p, bucketSize, method);
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

  init(lpMethod, false);
}
// INT DIRECT
Lpm::Lpm(
  const int t_pn,
  double* xx,
  const int t_N,
  const int p,
  const LpmMethod lpMethod,
  const int bucketSize,
  const int method
) {
  set_indirect = false;

  N = t_N;

  sample = new int[N];
  iprobabilities = new int[N];

  tree = new KDTree(xx, N, p, bucketSize, method);
  tree->init();

  if (N == 0 || t_pn == 0) {
    idx = new IndexList(0);
  } else if (t_pn == N) {
    idx = new IndexList(0);
    for (int i = 0; i < N; i++) {
      addUnitToSample(i);
    }
  } else {
    idx = new IndexList(N);
    for (int i = 0; i < N; i++) {
      probabilities[i] = t_pn;
      idx->set(i);
    }
  }

  init(lpMethod, true);
}
// DOUBLE INDIRECT
Lpm::Lpm(
  double* t_probabilities,
  KDTree* t_tree,
  IndexList* t_idx,
  const int t_N,
  const LpmMethod lpMethod,
  const double t_eps
) {
  set_indirect = true;
  sample = new int[N];

  probabilities = t_probabilities;
  tree = t_tree;
  idx = t_idx;
  N = t_N;
  eps = t_eps;

  init(lpMethod, false);
}
// INT INDIRECT
Lpm::Lpm(
  int* t_probabilities,
  KDTree* t_tree,
  IndexList* t_idx,
  const int t_N,
  const LpmMethod lpMethod
) {
  set_indirect = true;

  iprobabilities = t_probabilities;
  tree = t_tree;
  idx = t_idx;
  N = t_N;

  sample = new int[N];

  init(lpMethod, true);
}

void Lpm::init(const LpmMethod lpMethod, const bool isInt) {
  neighbours = new int[N];

  if (isInt) {
    _run = &Lpm::run_int;
  } else {
    _run = &Lpm::run_double;
  }

  set_run = true;

  switch (lpMethod) {
  case LpmMethod::lpm1:
    neighbours2 = new int[N];
    _draw = &Lpm::draw_lpm1;
    break;

  case LpmMethod::lpm2:
    _draw = &Lpm::draw_lpm2;
    break;

  case LpmMethod::lpm1search:
    neighbours2 = new int[N];
    history = new int[N];
    _draw = &Lpm::draw_lpm1search;
    break;

  case LpmMethod::rpm:
    _draw = &Lpm::draw_rpm;
    break;

  case LpmMethod::spm:
    _draw = &Lpm::draw_spm;
    break;

  default:
    std::invalid_argument("lpMethod does not exist");
    break;
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
  if (tree != nullptr)
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

void Lpm::draw_rpm(int* pair) {
  pair[0] = idx->draw();
  int len = idx->length() - 1;
  pair[1] = idx->get(intuniform(len));
  if (pair[0] == pair[1])
    pair[1] = idx->get(len);
  return;
}

void Lpm::draw_spm(int* pair) {
  if (!idx->exists(pair[0])) {
    pair[0] = pair[1];

    if (!idx->exists(pair[0])) {
      pair[0] += 1;

      if (pair[0] >= N)
        std::range_error("invalid value of pair 0");
    }
  }

  if (pair[0] >= pair[1]) {
    pair[1] = pair[0] + 1;
    return;
  }

  if (!idx->exists(pair[1])) {
    pair[1] += 1;
    if (pair[1] >= N)
      std::range_error("invalid value of pair 1");
  }

  return;
}

void Lpm::run_double() {
  int *pair = new int[2];
  pair[0] = 0;
  pair[1] = 1;

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

      if (pbig(*p1, eps))
        addUnitToSample(id1);
    }

    if (pclose(*p2, eps)) {
      eraseUnit(id2);

      if (pbig(*p2, eps))
        addUnitToSample(id2);
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);

    if (stduniform() < probabilities[id1])
      addUnitToSample(id1);

    eraseUnit(id1);
  }

  delete[] pair;

  std::sort(sample, sample + sampleSize);
  return;
}

void Lpm::run_int() {
  int *pair = new int[2];
  pair[0] = 0;
  pair[1] = 1;

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

      if (*p1 == N)
        addUnitToSample(id1);
    }

    if (*p2 == 0 || *p2 == N) {
      eraseUnit(id2);

      if (*p2 == N)
        addUnitToSample(id2);
    }
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);

    if (intuniform(N) < iprobabilities[id1])
      addUnitToSample(id1);

    eraseUnit(id1);
  }

  delete[] pair;

  std::sort(sample, sample + sampleSize);
  return;
}

