#include <algorithm>
// #include <Rcpp.h>

#include "CubeStratifiedClass.h"
#include "utils.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

CubeStratified::CubeStratified(
  double* t_prob,
  double* t_xbalance,
  int* t_strata,
  const int t_N,
  const int t_pb,
  const double t_eps
) {
  cubeMethod = CubeMethod::cube;

  Init(
    t_N,
    t_pb,
    t_eps,
    t_prob,
    t_xbalance,
    t_strata
  );
}

CubeStratified::CubeStratified(
  double* t_prob,
  double* t_xbalance,
  double* t_xspread,
  int* t_strata,
  const int t_N,
  const int t_pb,
  const int t_ps,
  const double t_eps,
  const int t_bucketSize,
  const int t_method
) {
  cubeMethod = CubeMethod::lcube;
  ps = t_ps;
  r_xspread = t_xspread;
  treeBucketSize = t_bucketSize;
  treeMethod = t_method;

  tree = new KDTree(r_xspread, t_N, ps, treeBucketSize, treeMethod);
  tree->init();
  spread = new double[t_N * ps];

  Init(
    t_N,
    t_pb,
    t_eps,
    t_prob,
    t_xbalance,
    t_strata
  );
}

void CubeStratified::Init(
  const int t_N,
  const int t_pb,
  const double t_eps,
  double* t_prob,
  double* t_xbalance,
  int* t_strata
) {
  N = t_N;
  pb = t_pb;
  eps = t_eps;
  r_prob = t_prob;
  r_xbalance = t_xbalance;
  r_strata = t_strata;

  idx = new IndexList(N);

  sample = std::vector<int>(0);
  sample.reserve(N);

  for (int i = N; i-- > 0; ) {
    idx->set(i);

    if (ProbabilityInt(r_prob[i], eps)) {
      EraseUnit(i);

      if (Probability1(r_prob[i], eps))
        AddUnitToSample(i);

      continue;
    }

    if (stratumMap.count(r_strata[i]) == 0)
      stratumMap[r_strata[i]] = 1;
    else
      stratumMap[r_strata[i]] += 1;
  }

  cube = new Cube(cubeMethod, N, pb + stratumMap.size(), eps);
  cube->InitIndirect();

  probabilities = std::vector<double>(N);
  index = std::vector<int>(0); index.reserve(N);
  stratumArr = std::vector<int>(stratumMap.size());
}

void CubeStratified::AddUnitToSample(const int id) {
  sample.push_back(id + 1);
  return;
}

void CubeStratified::EraseUnit(const int id) {
  idx->erase(id);
  if (tree != nullptr)
    tree->removeUnit(id);
  return;
}

void CubeStratified::Run() {
  RunFlightPerStratum();
  RunFlightOnFull();
  RunLandingPerStratum();

  std::sort(sample.begin(), sample.end());

  return;
}

void CubeStratified::RunFlightPerStratum() {
  int maxSize = pb + 2;
  cube->p = maxSize - 1;

  for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    // We can skip this completely if there are too few obs in the stratum
    if (it->second < maxSize)
      continue;

    cube->N = it->second;
    cube->sample.resize(0);
    cube->idx->resize(it->second);
    index.resize(0);

    // Add units into cube
    for (int i = 0; i < idx->length(); i++) {
      int id = idx->get(i);
      if (it->first != r_strata[id])
        continue;

      int indexSize = index.size();
      index.push_back(id);
      cube->idx->set(indexSize);
      cube->probabilities[indexSize] = r_prob[id];
      cube->amat[MatrixIdxRow(0, indexSize, it->second)] = 1.0;

      for (int k = 0; k < pb; k++)
        cube->amat[MatrixIdxRow(k + 1, indexSize, it->second)] =
          r_xbalance[MatrixIdxCol(id, k, N)] / r_prob[id];

      if (tree != nullptr) {
        for (int k = 0; k < ps; k++)
          spread[MatrixIdxRow(indexSize, k, ps)] =
            r_xspread[MatrixIdxRow(id, k, ps)];
      }
    }

    if (tree != nullptr) {
      KDTree* ttree = new KDTree(spread, it->second, ps, treeBucketSize, treeMethod);
      ttree->init();
      cube->tree = ttree;

      cube->RunFlight();

      delete ttree;
    } else {
      cube->RunFlight();
    }

    // Update idx and probabilities
    for (int i = 0; i < index.size(); i++) {
      int id = index[i];

      if (!cube->idx->exists(i)) {
        EraseUnit(id);
        it->second -= 1;
      } else {
        probabilities[id] = cube->probabilities[i];
      }
    }

    // Update sample
    for (int i = 0; i < cube->sample.size(); i++)
      AddUnitToSample(index[cube->sample[i] - 1]);

    // Add empty stratas for removal
    if (it->second == 0)
      stratumArr.push_back(it->first);
  }

  for (int k = 0; k < stratumArr.size(); k++)
    stratumMap.erase(stratumArr[k]);

  // Repopulate stratumArr with stratum numbers
  stratumArr.resize(0);
  for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it)
    stratumArr.push_back(it->first);

  return;
}

void CubeStratified::RunFlightOnFull() {
  int maxSize = pb + 1 + stratumMap.size();
  cube->p = maxSize - 1;

  if (idx->length() >= maxSize) {
    IndexList* tidx = cube->idx;
    cube->idx = idx;
    cube->tree = tree;
    int idxlen = idx->length();
    cube->N = N;
    cube->sample.resize(0);
    index.resize(0);

    for (int i = 0; i < idxlen; i++) {
      int id = idx->get(i);
      index.push_back(id);
      cube->probabilities[id] = probabilities[id];

      for (int k = 0; k < stratumMap.size(); k++)
        cube->amat[MatrixIdxRow(k, id, N)] = r_strata[id] == stratumArr[k] ? 1.0 : 0.0;

      for (int k = 0; k < pb; k++)
        cube->amat[MatrixIdxRow(stratumMap.size() + k, id, N)] =
          r_xbalance[MatrixIdxCol(id, k, N)] / r_prob[id];
          // r_xbalance[MatrixIdx(id, k, pb)] / r_prob[id];
    }

    cube->RunFlight();

    // Update sample
    for (int i = 0; i < cube->sample.size(); i++)
      AddUnitToSample(cube->sample[i] - 1);

    // Remove units from strata that has been decided
    for (int i = 0; i < idxlen; i++) {
      int id = index[i];
      if (!idx->exists(id))
        stratumMap[r_strata[id]] -= 1;
      else
        probabilities[id] = cube->probabilities[id];
    }

    // Put back tidx
    cube->idx = tidx;
  }

  return;
}

void CubeStratified::RunLandingPerStratum() {
  int maxSize = pb + 2;
  cube->p = maxSize - 1;
  cube->tree = nullptr;

  for (std::unordered_map<int, int>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    // Skip if all decided
    if (it->second <= 0)
      continue;

    cube->N = it->second;
    cube->sample.resize(0);
    cube->idx->resize(it->second);
    index.resize(0);

    // Add units into cube
    for (int i = 0; i < idx->length(); i++) {
      int id = idx->get(i);
      if (it->first != r_strata[id])
        continue;

      int indexSize = index.size();
      index.push_back(id);
      cube->idx->set(indexSize);
      cube->probabilities[indexSize] = probabilities[id];
      cube->amat[MatrixIdxRow(0, indexSize, it->second)] = 1.0;

      for (int k = 0; k < pb; k++)
        cube->amat[MatrixIdxRow(k + 1, indexSize, it->second)] =
          r_xbalance[MatrixIdxCol(id, k, N)] / probabilities[id];
          // r_xbalance[MatrixIdx(id, k, pb)] / probabilities[id];
    }

    cube->RunLanding();

    // Update sample
    for (int i = 0; i < cube->sample.size(); i++)
      AddUnitToSample(index[cube->sample[i] - 1]);
  }

  return;
}
