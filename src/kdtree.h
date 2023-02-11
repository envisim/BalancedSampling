#ifndef KDTREE_HEADER
#define KDTREE_HEADER

#include <float.h>
#include <algorithm>
#include <iostream>

class KDNode {
  // REGULAR NODE
public:
  KDNode *parent;
  KDNode *cleft; // <=
  KDNode *cright; // >
  int split;
  double value;
  // double *min;
  // double *max;

private:
  // TERMINAL NODE
  int terminal; // 1 if terminal node, 0 if regular node
public:
  int *units;
  int nunits;

public:
  KDNode(KDNode*, int, int);
  void setTerminal(int);
  int isTerminal();
  KDNode* getSibling();
  void reserveUnits(int);
  void addUnit(int);
  void removeUnit(int);
  int* getUnits() {return units;};
  int getNUnits() {return nunits;};
};

class KDTree {
private:
  double *data; // data array of length Np
  int N, p; // Data dimensions
  int bucketSize;

public:
  KDNode *top;

private:
  void splitSpread(KDNode*, int*, int, int);
  void splitMidpointSlide(KDNode*, int*, int, int);
  void findNeighbourInNode(int*, const int, int*, KDNode*, double*, int);
  void findNeighbourSearch(int*, const int, int*, KDNode*, double*, int);

public:
  KDTree(double*, int, int, int);
  void init();
  KDNode* findNode(int);
  int findNeighbour(int*, const int, const int);
  void removeUnit(int);
  double distance(double*, double*);
  double distanceIdx(int, int);
};

#endif
