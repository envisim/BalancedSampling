#ifndef LPM2INTERNAL_HEADER
#define LPM2INTERNAL_HEADER

#include "kdtree.h"
#include "uniform.h"
#include "index-list.h"

void lpm2_internal(
  KDTree *, //tree
  IndexList *, //idx
  double *, //probabilities
  int *, //neighbours
  const int, //N
  const double, // eps
  int *, //sample
  int * //sampleSize
);

#endif
