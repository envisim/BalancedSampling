#ifndef CPSINTERNAL_HEADER
#define CPSINTERNAL_HEADER

#include <algorithm>
#include <stdexcept>
#include "kdtree-cps.h"
#include "uniform.h"
#include "index-list.h"

enum CpsMethod {
  err = 0,
  lcps = 1,
  scps = 2,
  scpscoord = 3
};

CpsMethod intToCpsMethod(const int i);

class Cps {
protected:
  bool set_indirect = false;
  bool set_draw = false;

  int (Cps::*_draw)() = nullptr;
  double (Cps::*_random)(const int) = nullptr;

public:
  int N;
  double eps = 1e-12;

  IndexList* idx = nullptr;
  KDTreeCps* tree = nullptr;

  double* probabilities = nullptr;

  int* neighbours = nullptr;
  int* neighbours2 = nullptr;
  double* distances = nullptr;
  double* weights = nullptr;

  double* randomValues = nullptr;

  int histn = 0;

  int* sample = nullptr;
  int sampleSize = 0;

  // DIRECT
  Cps(
    const double* t_probabilities,
    double* xx,
    double* t_random,
    const int t_N,
    const int p,
    const CpsMethod cpsMethod,
    const int bucketSize,
    const int method,
    const double t_eps
  );
  // INDIRECT
  Cps(
    double* t_probabilities,
    KDTreeCps* t_tree,
    IndexList* t_idx,
    double* t_random,
    const int t_N,
    const CpsMethod cpsMethod,
    const double t_eps
  );

  ~Cps() {
    if (!set_indirect) {
      delete idx;
      delete tree;
      delete[] probabilities;
    }

    delete[] neighbours;
    delete[] neighbours2;
    delete[] distances;
    delete[] weights;
    delete[] sample;
  };

public:
  void init(const CpsMethod cpMethod);
  void addUnitToSample(const int id);
  void eraseUnit(const int id);
  void decideUnit(const int id);

protected:
  int draw_lcps();
  int draw_scps();
  int draw_scpscoord();

  double random_std(const int id);
  double random_arr(const int id);

  int draw() {return (this->*_draw)();};
  double random(const int id) {return (this->*_random)(id);};

  void _run();

public:
  void run() {
    if (!set_draw)
      std::runtime_error("_draw is nullptr");

    _run();
  };
};

#endif
