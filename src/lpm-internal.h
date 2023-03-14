#ifndef LPMINTERNAL_HEADER
#define LPMINTERNAL_HEADER

#include <functional>
#include "kdtree.h"
#include "uniform.h"
#include "index-list.h"

struct Lpm2Search {
  KDTree* tree;
  IndexList* idx;
  int* neighbours;
  int N;
  Lpm2Search(KDTree* t_tree, IndexList* t_idx, const int t_N) {
    tree = t_tree;
    idx = t_idx;
    neighbours = new int[t_N];
    N = t_N;
  };
  ~Lpm2Search() {
    delete[] neighbours;
  };
  void search(int* pair) {
    pair[0] = idx->draw();
    int len = tree->findNeighbour(neighbours, N, pair[0]);
    pair[1] = neighbours[intuniform(len)];
    return;
  };
};

void lpm_internal(
  KDTree*, //tree
  IndexList*, //idx
  double*, //probabilities
  const int, //N
  const double, // eps
  int*, //sample
  int*, //sampleSize
  std::function<void (int*)>
);

void lpm_int_internal(
  KDTree*, //tree
  IndexList*, //idx
  int*, //probabilities
  const int, //N
  int*, //sample
  int*, //sampleSize
  std::function<void (int*)>
);

#endif
