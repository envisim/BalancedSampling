#ifndef CUBEINTERNAL_HEADER
#define CUBEINTERNAL_HEADER

#include <cmath>
#include <stdexcept>
#include "kdtree.h"
#include "index-list.h"
#include "uniform.h"
#include "rref.h"

void cubeUpdatePhase(
  double*,
  double*,
  double*,
  int*,
  IndexList*,
  int*,
  int*,
  const double,
  const int
);

void cubeFlightPhase(
  const double*,
  double*,
  IndexList*,
  KDTree*,
  int*,
  int*,
  const double,
  const int,
  const int,
  int*,
  double*,
  double*,
  double*
);

void cubeLandingPhase(
  const double*,
  double*,
  IndexList*,
  /* KDTree*, */
  int*,
  int*,
  const double,
  const int,
  const int,
  int*,
  double*,
  double*
);

void cubeRunPhases(
  const double*,
  double*,
  IndexList*,
  KDTree*,
  int*,
  int*,
  const double,
  const int,
  const int
);

#endif
