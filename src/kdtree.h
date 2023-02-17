#ifndef KDTREE_HEADER
#define KDTREE_HEADER

#include <float.h>
#include <algorithm>
// #include <iostream>

class KDNode {
  // REGULAR NODE
public:
  KDNode *parent = nullptr;
  KDNode *cleft = nullptr; // <=
  KDNode *cright = nullptr; // >
  int split = -1;
  double value = 0.0;
  // double *min;
  // double *max;

private:
  // TERMINAL NODE
  int terminal = 0; // 1 if terminal node, 0 if regular node
  int nunits = 0;
public:
  int *units = nullptr;

public:
  KDNode(KDNode*, const int);
  ~KDNode();
  void copy(KDNode *original);
  void setTerminal(const int);
  int isTerminal();
  KDNode* getSibling();
  void addUnits(const int*, const int);
  void removeUnit(const int);
  int getNUnits();
  /* int* getUnits() {return units;}; */
};

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
  KDNode* findNode(const int);
  int findNeighbour(int*, const int, const int);
  int findClosest(int*, const int, const double*);
  void removeUnit(const int);
  double distanceIdx(const int, const int);
};

#endif
