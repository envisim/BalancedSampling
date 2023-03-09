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
