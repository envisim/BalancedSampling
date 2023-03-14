#ifndef SCPSINTERNAL_HEADER
#define SCPSINTERNAL_HEADER

#include <functional>
#include "kdtree-cps.h"
#include "index-list.h"

struct ScpsDecideProps {
  IndexList *idx;
  KDTreeCps *tree;
  double *probs;
  int *sample;
  int *sampleSize;
  double eps;

  ScpsDecideProps(
    IndexList *t_idx,
    KDTreeCps *t_tree,
    double *t_probs,
    int *t_sample,
    int *t_sampleSize,
    double t_eps
  ) :
  idx(t_idx),
    tree(t_tree),
    probs(t_probs),
    sample(t_sample),
    sampleSize(t_sampleSize),
    eps(t_eps)
    {};
};

void scpsDecide(
  ScpsDecideProps*,
  const int
);

void scps_internal(
  const int,
  double *,
  KDTreeCps *,
  const double,
  std::function<double (const int)>,
  std::function<int ()>,
  int *,
  int *,
  IndexList *
);

#endif
