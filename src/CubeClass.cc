#include <cmath>
#include <stdexcept>
#include <Rcpp.h>

#include "CubeClass.h"
#include "ReducedRowEchelonForm.h"
#include "uniform.h"
#include "utils.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

CubeMethod IntToCubeMethod(const int i) {
  if (1 <= i && i <= 2)
    return static_cast<CubeMethod>(i);

  std::invalid_argument("cube-method does not exist");
  return CubeMethod::err;
}

void Cube::InitIndirect() {
  sample = std::vector<int>(0);
  sample.reserve(N);
  probabilities = std::vector<double>(N);
  amat = std::vector<double>(N * p);

  index = std::vector<int>(N);
  uvec = std::vector<double>(p + 1);
  bmat = std::vector<double>((p + 1) * p);

  idx = new IndexList(N);

  switch(cubeMethod) {
  case CubeMethod::cube:
    _Draw = &Cube::Draw_cube;
    break;
  case CubeMethod::lcube:
    _Draw = &Cube::Draw_lcube;
    neighbours = std::vector<int>(N);
    distances = std::vector<double>(N);
    break;
  default:
    std::invalid_argument("cubeMethod does not exist");
    break;
  }

  set_draw = true;
};

void Cube::Init(
  const double* t_probabilities,
  double* xxbalance,
  const int t_N,
  const int t_p,
  const double t_eps
) {
  N = t_N;
  p = t_p;
  eps = t_eps;

  if (N == 0) {
    idx = new IndexList(0);
    return;
  }

  InitIndirect();

  // Decrement done before evaluating the loop, so it begins on N - 1
  for (int i = N; i-- > 0; ) {
    probabilities[i] = t_probabilities[i];
    idx->set(i);

    if (ProbabilityInt(probabilities[i], eps)) {
      EraseUnit(i);

      if (Probability1(probabilities[i], eps))
        AddUnitToSample(i);

      continue;
    }

    for (int k = 0; k < p; k++)
      amat[MatrixIdxRow(k, i, N)] =
        xxbalance[MatrixIdxCol(i, k, p)] / probabilities[i];
  }

  return;
};

// NO INIT CUBE
Cube::Cube(
  const CubeMethod t_cubeMethod,
  const int t_N,
  const int t_p,
  const double t_eps
) {
  cubeMethod = t_cubeMethod;
  N = t_N;
  p = t_p;
  eps = t_eps;
  set_indirect = true;
}

// DIRECT CUBE
Cube::Cube(
  const double* t_probabilities,
  double* xxbalance,
  const int t_N,
  const int t_p,
  const double t_eps
) {
  cubeMethod = CubeMethod::cube;
  set_indirect = false;
  Init(t_probabilities, xxbalance, t_N, t_p, t_eps);

  idx->shuffle();
}

// DIRECT LCUBE
Cube::Cube(
  const double* t_probabilities,
  double* xxbalance,
  const int t_N,
  const int t_pbalance,
  const double t_eps,
  double* xxspread,
  const int t_pspread,
  const int bucketSize,
  const int method
) {
  cubeMethod = CubeMethod::lcube;
  set_indirect = false;
  Init(t_probabilities, xxbalance, t_N, t_pbalance, t_eps);

  tree = new KDTree(xxspread, N, t_pspread, bucketSize, method);
  tree->init();
}

void Cube::AddUnitToSample(const int id) {
  sample.push_back(id + 1);
  return;
}

void Cube::EraseUnit(const int id) {
  idx->erase(id);

  // Needed like this as tree might be nullptr during landing
  if (tree != nullptr)
    tree->removeUnit(id);

  return;
}

int Cube::MaxSize() {
  int il = idx->length();
  return p + 1 <= il ? p + 1 : il;
}

void Cube::Draw(const int maxSize) {return (this->*_Draw)(maxSize);}

void Cube::RunFlight() {
  if (!set_draw)
    std::runtime_error("_Draw is nullptr");

  _RunFlight();
}

void Cube::RunLanding() {
  if (!set_draw)
    std::runtime_error("_Draw is nullptr");

  _RunLanding();
}

void Cube::Run() {
  if (!set_draw)
    std::runtime_error("_Draw is nullptr");

  _RunFlight();
  _RunLanding();
  std::sort(sample.begin(), sample.end());

  return;
}

void Cube::RunUpdate() {
  int maxSize = MaxSize();
  ReducedRowEchelonForm(&bmat[0], maxSize - 1, maxSize);

  double lambda1 = DBL_MAX;
  double lambda2 = DBL_MAX;

  for (int i = 0; i < maxSize; i++) {
    if (i == maxSize - 1) {
      uvec[i] = 1.0;
    } else {
      uvec[i] = -bmat[MatrixIdxRow(i, maxSize - 1, maxSize)];
    }

    double lval1 = std::abs(probabilities[index[i]] / uvec[i]);
    double lval2 = std::abs((1.0 - probabilities[index[i]]) / uvec[i]);

    if (uvec[i] >= 0.0) {
      if (lambda1 > lval2)
        lambda1 = lval2;
      if (lambda2 > lval1)
        lambda2 = lval1;
    } else {
      if (lambda1 > lval1)
        lambda1 = lval1;
      if (lambda2 > lval2)
        lambda2 = lval2;
    }
  }

  double lambda = stduniform(lambda1 + lambda2) < lambda2 ? lambda1 : -lambda2;

  for (int i = 0; i < maxSize; i++) {
    int id = index[i];
    probabilities[id] += lambda * uvec[i];

    if (ProbabilityInt(probabilities[id], eps)) {
      EraseUnit(id);

      if (Probability1(probabilities[id], eps))
        AddUnitToSample(id);
    }
  }

  return;
}

void Cube::Draw_cube(const int maxSize) {
  for (int i = 0; i < maxSize; i++)
    index[i] = idx->get(i);

  return;
}

void Cube::Draw_lcube(const int maxSize) {
  int id = idx->draw();
  index[0] = id;
  int len = tree->findNeighboursN(&distances[0], &index[1], maxSize - 1, id);
  len += 1;

  if (len > maxSize) {
    // We have equals at the end
    int i = len - 1;
    while (i > 2) {
      if (distances[index[i]] > distances[index[i - 1]])
        break;

      i -= 1;
    }

    for (; i < maxSize; i++) {
      int j = intuniform(len - i) + i;
      if (i == j)
        continue;
      int temp = index[i];
      index[i] = index[j];
      index[j] = temp;
    }
  }

  return;
}

void Cube::_RunFlight() {
  // Cases:
  // - choose from tree and all
  // - choose from list and all

  int maxSize = MaxSize();

  while (idx->length() >= maxSize) {
    Draw(maxSize);

    // Prepare bmat
    for (int i = 0; i < maxSize; i++) {
      for (int k = 0; k < maxSize - 1; k++) {
        bmat[MatrixIdxRow(k, i, maxSize)] = amat[MatrixIdxRow(k, index[i], N)];
      }
    }

    RunUpdate();
  }

  return;
}

void Cube::_RunLanding() {
  if (idx->length() >= p + 1)
    std::range_error("landingphase committed early");

  while (idx->length() > 1) {
    int maxSize = idx->length();

    for (int i = 0; i < maxSize; i++) {
      index[i] = idx->get(i);

      for (int k = 0; k < maxSize - 1; k++) {
        bmat[MatrixIdxRow(k, i, maxSize)] = amat[MatrixIdxRow(k, index[i], N)];
      }
    }

    RunUpdate();
  }

  if (idx->length() == 1) {
    int id1 = idx->get(0);

    if (stduniform() < probabilities[id1])
      AddUnitToSample(id1);

    EraseUnit(id1);
  }

  return;
}
