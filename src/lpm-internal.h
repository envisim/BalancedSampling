#ifndef LPMINTERNAL_HEADER
#define LPMINTERNAL_HEADER

#include <algorithm>
#include <functional>
#include <stdexcept>
#include "kdtree.h"
#include "uniform.h"
#include "index-list.h"

enum LpmMethod {
  err = 0,
  lpm1 = 1,
  lpm2 = 2,
  lpm1search = 3,
  rpm = 4,
  spm = 5
};

LpmMethod intToLpmMethod(const int i);

class Lpm {
protected:
  bool set_indirect = false;
  bool set_draw = false;
  bool set_run = false;

  void (Lpm::*_draw)(int*) = nullptr;
  void (Lpm::*_run)() = nullptr;

public:
  int N;
  double eps = 1e-12;

  IndexList* idx = nullptr;
  KDTree* tree = nullptr;

  double* probabilities = nullptr;
  int* iprobabilities = nullptr;

  int* neighbours = nullptr;
  int* neighbours2 = nullptr;

  int* history = nullptr;
  int histn = 0;

  int* sample = nullptr;
  int sampleSize = 0;

// DOUBLE
  Lpm(
    const double* t_probabilitites,
    double* xx,
    const int t_N,
    const int p,
    const LpmMethod lpMethod,
    const int bucketSize,
    const int method,
    const double t_eps
  );
  // INT
  Lpm(
    const int t_pn,
    double* xx,
    const int t_N,
    const int p,
    const LpmMethod lpMethod,
    const int bucketSize,
    const int method
  );
  // DOUBLE PRE
  Lpm(
    double* t_probabilities,
    KDTree* t_tree,
    IndexList* t_idx,
    const int t_N,
    const LpmMethod lpMethod,
    const double t_eps
  );
  // INT PRE
  Lpm(
    int* t_probabilities,
    KDTree* t_tree,
    IndexList* t_idx,
    const int t_N,
    const LpmMethod lpMethod
  );

  ~Lpm() {
    if (!set_indirect) {
      delete idx;
      delete tree;
      delete[] probabilities;
      delete[] iprobabilities;
    }

    delete[] neighbours;
    delete[] neighbours2;
    delete[] history;
    delete[] sample;
  };

public:
  void init(const LpmMethod lpMethod, const bool isInt);
  void addUnitToSample(const int id);
  void eraseUnit(const int id);

protected:
  void draw_lpm1(int* pair);
  void draw_lpm2(int* pair);
  void draw_lpm1search(int* pair);
  void draw_rpm(int *pair);
  void draw_spm(int *pair);
  void run_double();
  void run_int();

  void draw(int* pair) {(this->*_draw)(pair);};

public:
  void run() {
    if (!set_run)
      std::runtime_error("_run is nullptr");
    if (!set_draw)
      std::runtime_error("_draw is nullptr");

    (this->*_run)();
  };
  // Return
};

#endif
