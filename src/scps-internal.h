#ifndef SCPSINTERNAL_HEADER
#define SCPSINTERNAL_HEADER

#include <functional>
#include "kdtree-cps.h"
#include "index-list.h"

struct ScpsDecideProps {
  KDTreeCps* tree;
  IndexList* idx;
  double* probs;
  double eps;
  int* sample;
  int* sampleSize;

  ScpsDecideProps(
    KDTreeCps* t_tree,
    IndexList* t_idx,
    double* t_probs,
    double t_eps,
    int* t_sample,
    int* t_sampleSize
  ) :
  tree(t_tree),
    idx(t_idx),
    probs(t_probs),
    eps(t_eps),
    sample(t_sample),
    sampleSize(t_sampleSize)
    {};
};

void scpsDecide(
  ScpsDecideProps*,
  const int
);

void scps_internal(
  KDTreeCps*,
  IndexList*,
  double*,
  const int,
  const double,
  int*,
  int*,
  std::function<double (const int)>,
  std::function<int ()>
);

#endif
