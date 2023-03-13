#ifndef KDTREECPS_HEADER
#define KDTREECPS_HEADER

#include <float.h>
#include <algorithm>
#include "kdtree.h"

class KDTreeCps : public KDTree {
public:
  using KDTree::KDTree;
  int findNeighbours(const double*, double*, double*, int*, const int);
protected:
  void findNeighboursSearch(const double*, double*, double*, int*, int*, KDNode*, double*, const int, double*);
  void findNeighboursInNode(const double*, double*, double*, int*, int*, KDNode*, double*, const int, double*);
};

#endif
