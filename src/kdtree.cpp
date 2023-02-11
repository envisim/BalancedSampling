#include "kdtree.h"

double square(double value) {
  return value * value;
}

KDNode::KDNode(KDNode *par, int term, int bs) {
  cleft = nullptr;
  cright = nullptr;
  units = nullptr;
  nunits = 0;

  parent = par;
  setTerminal(term);
  if (terminal == 1) {
    reserveUnits(bs);
  }

  // min = new double[p];
  // max = new double[p];
}

void KDNode::setTerminal(int term) {
  // terminal = term > 0 ? 1 : 0;
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

void KDNode::reserveUnits(int len) {
  delete[] units;
  units = new int[len];
  nunits = 0;
}

void KDNode::addUnit(int idx) {
  units[nunits] = idx;
  nunits += 1;
}

void KDNode::removeUnit(int idx) {
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

KDTree::KDTree(double *dt, int NN, int pp, int bs) {
  data = dt;
  N = NN;
  p = pp;
  bucketSize = bs;
  top = nullptr;
}

void KDTree::init() {
  if (bucketSize < 1 || N < 1 || p < 1)
    return;

  top = new KDNode(nullptr, bucketSize - N + 1, bucketSize);

  if (top->isTerminal() == 1) {
    for (int i = 0; i < N; i++)
      top->addUnit(i);
    return;
  }

  int *units = new int[N];
  for (int i = 0; i < N; i++)
    units[i] = i;

  splitSpread(top, units, 0, N);

  delete[] units;
}

void KDTree::splitSpread(KDNode *parent, int *units, int l, int r) {
  double *mins = new double[p];
  double *maxs = new double[p];

  std::copy_n(data + units[l] * p, p, mins);
  std::copy_n(data + units[l] * p, p, maxs);

  for (int i = l + 1; i < r; i++) {
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

  double minval = mins[parent->split];
  double maxval = maxs[parent->split];
  delete [] mins;
  delete[] maxs;

  // We must handle if there is not information left, even though
  // the bucketSize is too small
  if (spread == 0.0) {
    parent->setTerminal(1);
    parent->reserveUnits(r - l);
    for (int *i = units + l; i != units + r; i++)
      parent->addUnit(*i);
    return;
  }

  double *dt = data + parent->split;

  std::sort(
    units + l,
    units + r,
    [&, parent, dt](const int &a, const int &b) {
      return *(dt + a * p) < *(dt + b * p);
      // return data[a * p + parent->split] < data[b * p + parent->split];
    }
  );

  int m = ((r - l) >> 1) + l;
  parent->value = *(dt + units[m - 1] * p);

  // If we chose a bad midpoint, we need to look further
  if (parent->value == *(dt + units[m] * p)) {
    int searchAbove = parent->value != maxval;
    int searchBelow = parent->value != minval;

    for (int i = 1; i < m - l; i++) {
      // Check above
      if (searchAbove && parent->value != *(dt + units[m + i] * p)) {
        m += i;
        break;
      }

      // Check below
      // May cause oob, but should not as spread > 0.0 and above is checked first
      // if l = 0, r = 5, m = 2: 0 1 2 3 4
      // we check 3, 0, 4 (1, 2 is checked at if)
      // and one must be different as spread > 0.0
      // Also, if 0 == value, then we don't check below
      //
      // if l = 0, r = 6, m = 3: 0 1 2 3 4 5
      // we check 4, 1, 5, 0 (2, 3 is checked at if)
      // and one must be different as spread > 0.0
      if (searchBelow && parent->value != *(dt + units[m - 1 - i] * p)) {
        m -= i;
        break;
      }
    }

    parent->value = *(dt + units[m - 1] * p);
  }

  // parent->value = *(dt + units[m - 1] * p);
  parent->cleft = new KDNode(parent, bucketSize - (m - l) + 1, bucketSize);
  parent->cright = new KDNode(parent, bucketSize - (r - m) + 1, bucketSize);

  if (parent->cleft->isTerminal()) {
    for (int *i = units + l; i != units + m; i++)
      parent->cleft->addUnit(*i);
  } else {
    splitSpread(parent->cleft, units, l, m);
  }

  if (parent->cright->isTerminal()) {
    for (int *i = units + m; i != units + r; i++)
      parent->cright->addUnit(*i);
  } else {
    splitSpread(parent->cright, units, m, r);
  }
}

// void KDTree::splitMidpointSlide(KDNode *parent, int *units, int l, int r) {
//   return;
// }

KDNode* KDTree::findNode(int idx) {
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

  int size = 0;
  double mindist = DBL_MAX;

  findNeighbourSearch(neighbours, maxsize, &size, top, &mindist, idx);

  return size;
}

void KDTree::findNeighbourSearch(
    int *neighbours, const int maxsize, int *size, KDNode *node, double *mindist, int idx
) {
  if (node->isTerminal()) {
    findNeighbourInNode(neighbours, maxsize, size, node, mindist, idx);
    return;
  }

  KDNode *nextnode = data[idx * p + node->split] <= node->value
    ? node->cleft
    : node->cright;

  findNeighbourSearch(neighbours, maxsize, size, nextnode, mindist, idx);

  // <= to be precise, but maybe we can ignore this...
  if (square(data[idx * p + node->split] - node->value) < *mindist) {
    findNeighbourSearch(neighbours, maxsize, size, nextnode->getSibling(), mindist, idx);
  }
}


void KDTree::findNeighbourInNode(
  int *neighbours, const int maxsize, int *size, KDNode *node, double *mindist, int idx
) {
  for (int i = 0; i < node->nunits; i++) {
    if (node->units[i] == idx) // Same unit
      continue;

    double dist = distanceIdx(idx, node->units[i]);

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

void KDTree::removeUnit(int idx) {
  KDNode *node = findNode(idx);
  node->removeUnit(idx);
  return;
}

double KDTree::distance(double *u1, double *u2) {
  double dist = 0.0;

  for (int k = 0; k < p; k++) {
    double temp = *(u1 + k) - *(u2 + k);
    dist += temp * temp;
  }

  return dist;
}

double KDTree::distanceIdx(int id1, int id2) {
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
