// #include <Rcpp.h>
#include "kdtree.h"

// KDNode::KDNode(KDNode *par, int term, int bs) {
KDNode::KDNode(KDNode *par, const int term) {
  parent = par;
  setTerminal(term);

  // min = new double[p];
  // max = new double[p];
}

KDNode::~KDNode() {
  delete[] units;
  delete cleft;
  delete cright;
  // delete this;
}

void KDNode::copy(KDNode *original) {
  if (terminal) {
    addUnits(original->units, original->nunits);
    return;
  }

  split = original->split;
  value = original->value;

  cleft = new KDNode(this, original->cleft->isTerminal());
  cleft->copy(original->cleft);

  cright = new KDNode(this, original->cright->isTerminal());
  cright->copy(original->cright);
}

void KDNode::setTerminal(const int term) {
  terminal = term > 0;
}

int KDNode::isTerminal() {
  return terminal;
}

KDNode* KDNode::getSibling() {
  if (parent == nullptr)
    return nullptr;

  return this == parent->cleft ? parent->cright : parent->cleft;
}

int KDNode::getNUnits() {
  return nunits;
}

void KDNode::addUnits(const int *aunits, const int n) {
  delete[] units;
  units = new int[n];

  for (int i = 0; i < n; i++)
    units[i] = aunits[i];

  nunits = n;
}

void KDNode::removeUnit(const int idx) {
  if (nunits == 0)
    return;

  int *it = std::find(units, units + nunits, idx);

  if (it == units + nunits)
    return;

  nunits -= 1;

  if (nunits > 0) {
    if (it != units + nunits)
      *it = units[nunits];
    return;
  }

  // ======== IF WE REMOVE TREES
  // KDNode *sibling = getSibling();

  // if (parent == nullptr || sibling == nullptr)
  //   return;

  // if (sibling->isTerminal()) {
  //   parent->setTerminal(1);
  //   parent->units = sibling->units;
  //   parent->nunits = sibling->nunits;
  // } else {
  //   parent->cleft = sibling->cleft;
  //   parent->cright = sibling->cright;
  //   parent->split = sibling->split;
  //   parent->value = sibling->value;
  //   parent->cleft->parent = parent;
  //   parent->cright->parent = parent;
  // }

  // delete sibling;
  // delete[] units;
  // delete this;
}

// ############################# K-D-TREE ######################################

KDTree::KDTree(double *dt, const int NN, const int pp, const int bs, const int sm) {
  data = dt;
  N = NN;
  p = pp;
  bucketSize = bs;
  top = nullptr;
  liml = new double[pp];
  limr = new double[pp];

  for (int i = 0; i < pp; i++) {
    liml[i] = DBL_MAX;
    limr[i] = -DBL_MAX;
  }

  if (sm == 0) {
    method = 0;
    splitMethod = &KDTree::splitMod;
  } else if (sm == 1) {
    method = 1;
    splitMethod = &KDTree::splitSpread;
  } else {
    method = 2;
    splitMethod = &KDTree::splitMidpointSlide;
  }
}

KDTree::~KDTree() {
  delete top;
  delete[] liml;
  delete[] limr;
  // delete this;
}

KDTree* KDTree::copy() {
  KDTree* nt = new KDTree(data, N, p, bucketSize, method);
  nt->top = new KDNode(nullptr, top->isTerminal());
  nt->top->copy(top);
  return nt;
}

void KDTree::init() {
  if (bucketSize < 1 || N < 1 || p < 1)
    return;

  top = new KDNode(nullptr, bucketSize - N + 1);

  int *units = new int[N];
  double *dt = data;
  for (int i = 0; i < N; i++) {
    units[i] = i;

    for (int j = 0; j < p; j++) {
      if (*dt < liml[j])
        liml[j] = *dt;
      if (*dt > limr[j])
        limr[j] = *dt;
      dt++;
    }
  }

  if (top->isTerminal() == 1) {
    top->addUnits(units, N);
    delete[] units;
    return;
  }

  split(top, units, N);

  delete[] units;
}

void KDTree::split(KDNode *parent, int *units, const int n) {
  // int m = splitMod(parent, units, n);
  // int m = splitSpread(parent, units, n);
  // int m = splitMidpointSlide(parent, units, n);
  int m = (this->*splitMethod)(parent, units, n);

  // Check if m == -1, then it is not possible to split any more, and
  // we need to accept the current bucket, no matter the size.
  if (m == -1) {
    parent->setTerminal(1);
    parent->addUnits(units, n);
    return;
  }

  parent->cleft = new KDNode(parent, bucketSize - m + 1);
  parent->cright = new KDNode(parent, bucketSize - (n - m) + 1);

  if (parent->cleft->isTerminal()) {
    parent->cleft->addUnits(units, m);
  } else {
    split(parent->cleft, units, m);
  }

  if (parent->cright->isTerminal()) {
    parent->cright->addUnits(units + m, n - m);
  } else {
    split(parent->cright, units + m, n - m);
  }
}

int KDTree::splitM(int *units, const int n, const int m0, const int k) {
  int *tunits = new int[n];
  double *dt = data + k;
  int l = 0;
  int r = n;
  double value = *(dt + units[m0] * p);

  for (int i = 0; i < n; i++) {
    double temp = *(dt + units[i] * p);
    if (temp < value) {
      tunits[l] = units[i];
      l += 1;
    } else if (temp > value) {
      r -= 1;
      tunits[r] = units[i];
    }
  }

  int m = l;

  for (int i = 0; i < n; i++) {
    double temp = *(dt + units[i] * p);
    if (temp == value) {
      tunits[m] = units[i];
      m += 1;
    }
  }

  for (int i = 0; i < n; i++) {
    units[i] = tunits[i];
  }

  delete[] tunits;

  if (m0 < l) {
    return splitM(units, l, m0, k);
  } else if (m0 >= r) {
    return r + splitM(units + r, n - r, m0 - r, k);
  }

  return m;
}

int KDTree::splitSpread(KDNode *parent, int *units, const int n) {
  double *mins = new double[p];
  double *maxs = new double[p];

  std::copy_n(data + units[0] * p, p, mins);
  std::copy_n(data + units[0] * p, p, maxs);

  for (int i = 0; i < n; i++) {
    for (int k = 0; k < p; k++) {
      double *temp = data + units[i] * p + k;
      if (*temp < mins[k]) {
        mins[k] = *temp;
      } else if (*temp > maxs[k]) {
        maxs[k] = *temp;
      }
    }
  }

  parent->split = 0;
  double spread = maxs[0] - mins[0];
  for (int k = 1; k < p; k++) {
    double temp = maxs[k] - mins[k];
    if (temp > spread) {
      parent->split = k;
      spread = temp;
    }
  }

  delete[] mins;
  delete[] maxs;

  if (spread == 0.0)
    return -1;

  int m = splitM(units, n, n >> 1, parent->split);

  parent->value = *(data + units[m - 1] * p + parent->split);

  return m;
}

int KDTree::splitMod(KDNode *parent, int *units, const int n) {
  int lvl = 0;
  KDNode *par = parent;
  while (par->parent != nullptr) {
    par = par->parent;
    lvl++;
  }

  parent->split = lvl % p;

  int m = splitM(units, n, n >> 1, parent->split);

  parent->value = *(data + units[m - 1] * p + parent->split);

  return m;
}

int KDTree::splitMidpointSlide(KDNode *parent, int *units, const int n) {
  double *mins = new double[p];
  double *maxs = new double[p];

  std::copy_n(liml, p, mins);
  std::copy_n(limr, p, maxs);

  KDNode *par = parent;
  while(par->parent != nullptr) {
    if (par->parent->cleft == par) {
      if (par->parent->value < maxs[par->parent->split])
        maxs[par->parent->split] = par->parent->value;
    } else {
      if (par->parent->value > mins[par->parent->split])
        mins[par->parent->split] = par->parent->value;
    }

    par = par->parent;
  }

  parent->split = 0;
  double spread = maxs[0] - mins[0];
  for (int k = 1; k < p; k++) {
    double temp = maxs[k] - mins[k];
    if (temp > spread) {
      parent->split = k;
      spread = temp;
    }
  }

  if (spread == 0.0)
    return -1;

  parent->value = (maxs[parent->split] + mins[parent->split]) * 0.5;

  delete[] mins;
  delete[] maxs;

  double *dt = data + parent->split;
  int l = 0;
  int r = n;
  double small = DBL_MAX;
  double big = -DBL_MAX;

  int *tunits = new int[n];

  for (int i = 0; i < n; i++) {
    double temp = *(dt + units[i] * p);
    if (temp <= parent->value) {
      tunits[l] = units[i];
      l += 1;

      if (temp > big)
        big = temp;
    } else {
      r -= 1;
      tunits[r] = units[i];

      if (temp < small)
        small = temp;
    }
  }

  if (l > 0 && r < n) {
    for (int i = 0; i < n; i++)
      units[i] = tunits[i];

    delete[] tunits;
    return l;
  }

  delete[] tunits;

  if (l == 0) {
    // If there are no units lesseq than value, we
    // FIND ALL SMALLEST
    for (int i = 0; i < n; i++) {
      double temp = *(dt + units[i] * p);
      if (temp == small) {
        if (i != l) {
          int t = units[l];
          units[l] = units[i];
          units[i] = t;
        }
        l += 1;
      }
    }

    return l;
  }

  if (r == n) {
    // If there are no units bigger than value, we
    // FIND ALL BIGGEST
    for (int i = n - 1; i >= 0; i--) {
      double temp = *(dt + units[i] * p);
      if (temp == big) {
        r -= 1;
        if (i != r) {
          int t = units[r];
          units[r] = units[i];
          units[i] = t;
        }
      }
    }

    return r;
  }


  return -1;
}

KDNode* KDTree::findNode(const int idx) {
  double *obs = data + idx * p;
  KDNode *node = top;

  while (node->isTerminal() == 0) {
    if (*(obs + node->split) <= node->value) {
      node = node->cleft;
    } else {
      node = node->cright;
    }
  }

  return node;
}

int KDTree::findNeighbour(int *neighbours, const int maxsize, const int idx) {
  if (top == nullptr)
    return 0;

  double *unit = data + idx * p;
  int size = 0;
  double mindist = DBL_MAX;

  findNeighbourSearch(neighbours, maxsize, &size, top, &mindist, idx, unit);

  return size;
}

void KDTree::findNeighbourSearch(
  int *neighbours, const int maxsize, int *size, KDNode *node, double *mindist, const int idx, const double *unit
) {
  if (node->isTerminal()) {
    findNeighbourInNode(neighbours, maxsize, size, node, mindist, idx, unit);
    return;
  }

  // KDNode *nextnode = data[idx * p + node->split] <= node->value
  //   ? node->cleft
  //   : node->cright;
  double dist = *(unit + node->split) - node->value;
  KDNode *nextnode = dist <= 0.0 ? node->cleft : node->cright;

  findNeighbourSearch(neighbours, maxsize, size, nextnode, mindist, idx, unit);

  // <= to be precise, but maybe we can ignore this...
  // if (square(data[idx * p + node->split] - node->value) < *mindist) {
  if (dist * dist < *mindist) {
    findNeighbourSearch(neighbours, maxsize, size, nextnode->getSibling(), mindist, idx, unit);
  }
}


void KDTree::findNeighbourInNode(
  int *neighbours, const int maxsize, int *size, KDNode *node, double *mindist, const int idx, const double *unit
) {
  int nunits = node->getNUnits();

  for (int i = 0; i < nunits; i++) {
    if (node->units[i] == idx) // Same unit
      continue;

    // double dist = distanceIdx(idx, node->units[i]);
    double dist = distance(unit, data + node->units[i] * p);

    if (dist < *mindist) {
      neighbours[0] = node->units[i];
      *size = 1;
      *mindist = dist;
    } else if (dist == *mindist && *size < maxsize) {
      neighbours[*size] = node->units[i];
      *size += 1;
    }
  }
}

void KDTree::removeUnit(const int idx) {
  KDNode *node = findNode(idx);
  node->removeUnit(idx);
  return;
}

double KDTree::distance(const double *u1, const double *u2) {
  double dist = 0.0;

  for (int k = 0; k < p; k++) {
    double temp = *(u1 + k) - *(u2 + k);
    dist += temp * temp;
  }

  return dist;
}

double KDTree::distanceIdx(const int id1, const int id2) {
  double *u1 = data + id1 * p;
  double *u2 = data + id2 * p;
  return distance(u1, u2);
}

#ifdef CLI

void printKDNode(KDNode *node, int par) {
  if (node->isTerminal() == 1) {
    for (int i = 0; i < par; i++)
      std::cout << "-";
    std::cout << " ";

    int nunits = node->getNUnits();
    int *units = node->getUnits();
    for (int i = 0; i < nunits; i++)
      std::cout << units[i] << " ";
    std::cout << std::endl;
    return;
  }

  for (int i = 0; i < par; i++)
    std::cout << "-";
  std::cout << " ";

  std::cout << "SPLIT at " << node->split << " (" << node->value << ")" << std::endl;
  printKDNode(node->cleft, par+1);
  printKDNode(node->cright, par+1);
}

void printKDTree(KDTree *tree) {
  printKDNode(tree->top, 0);
}

int main() {
  double x[20] = {
    0.68512126, 0.3251399, // 0 8
    0.05171296, 0.7740967, // 1 7
    0.74261974, 0.9242472, // 2 4
    0.13036790, 0.3870030, // 3 6
    0.77980495, 0.7918827, // 4 2
    0.01413735, 0.5849822, // 5 1
    0.25770368, 0.4773944, // 6 3
    0.09543018, 0.8095111, // 7 1
    0.39014922, 0.3908506, // 8 6
    0.32050716, 0.1994035  // 9 8
  };

  KDTree *tree = new KDTree(&x[0], 10, 2, 1);
  tree->init();
  printKDTree(tree);

  int *neighbours = new int[10];

  for (int i = 0; i < 10; i++) {
    int size = tree->findNeighbour(neighbours, 10, i);
    // tree->removeUnit(i);

    std::cout << i << ": ";

    for (int j = 0; j < size; j++)
      std::cout << neighbours[j] << " ";
    std::cout << std::endl;

    // printKDTree(tree);
  }

  return 0;
}

#endif
