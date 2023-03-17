#include <float.h>
#include <algorithm>
#include <stdexcept>
// #include <Rcpp.h>
#include "kdnode.h"

KDNode::KDNode(KDNode *par, const int term) {
  parent = par;
  setTerminal(term);
}

KDNode::~KDNode() {
  delete[] units;
  delete cleft;
  delete cright;
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

void KDNode::prune(const int bucketSize) {
  // We can't prune a terminal node
  if (terminal)
    return;

  // If a node is not terminal, but is missing a child, something went wrong
  if (cleft == nullptr || cright == nullptr)
    return;

  if (!cleft->terminal)
    cleft->prune(bucketSize);

  if (!cright->terminal)
    cleft->prune(bucketSize);

  if (!cleft->terminal || !cright->terminal)
    return;

  // Both nodes are terminal
  int cnunits = cleft->nunits + cright->nunits;
  if (cnunits > bucketSize)
    return;

  // The nodes are too small, prune
  delete[] units;
  units = new int[cnunits];
  nunits = 0;

  for (int i = 0; i < cleft->nunits && nunits < cnunits; i++) {
    units[nunits] = cleft->units[i];
    nunits += 1;
  }

  for (int i = 0; i < cright->nunits && nunits < cnunits; i++) {
    units[nunits] = cright->units[i];
    nunits += 1;
  }

  terminal = 1;
  delete cleft;
  delete cright;
  cleft = nullptr;
  cright = nullptr;
  return;
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

int KDNode::getSize() {
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
  // If there are no units here, nothing to remove
  if (nunits == 0)
    return;

  // Find the pointer to idx in units
  int *it = std::find(units, units + nunits, idx);

  // If nothing found, return
  if (it == units + nunits)
    return;

  // Something found, thus we have one less unit
  nunits -= 1;

  // If the iterator is pointing to the pos at nunits, we are done
  if (it == units + nunits)
    return;

  *it = units[nunits];
  return;
}

bool KDNode::exists(const int idx) {
  if (nunits == 0)
    return false;

  int *it = std::find(units, units + nunits, idx);
  return it != units + nunits;
}
