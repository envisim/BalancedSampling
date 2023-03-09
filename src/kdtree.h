#ifndef KDTREE_HEADER
#define KDTREE_HEADER

#include <float.h>
#include <algorithm>
#include <stdexcept>
#include "kdnode.h"

class KDTree {
protected:
  double *data; // data array of length Np
  int N, p; // Data dimensions
  int bucketSize;
  double *liml = nullptr;
  double *limr = nullptr;
  int method = 2;
  int (KDTree::*splitMethod)(KDNode*, int*, const int) = nullptr;

public:
  KDNode *top = nullptr;

protected:
  void split(KDNode*, int*, const int);
  int splitM(int*, const int, const int, const int);
  int splitSpread(KDNode*, int*, const int);
  int splitMod(KDNode*, int*, const int);
  int splitMidpointSlide(KDNode*, int*, const int);
  void findNeighbourInNode(int*, const int, int*, KDNode*, double*, const int, const double*);
  void findNeighbourSearch(int*, const int, int*, KDNode*, double*, const int, const double*);
  double distance(const double*, const double*);

public:
  KDTree(double*, const int, const int, const int, const int);
  ~KDTree();
  KDTree* copy();
  void init();
  void prune();

  KDNode* findNode(const int);
  int findNeighbour(int*, const int, const int);
  int findClosest(int*, const int, const double*);
  double findSmallestDistanceToPoint(const double*);
  void removeUnit(const int);
  bool unitExists(const int);
  double distanceIdx(const int, const int);
};

#endif
