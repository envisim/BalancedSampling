#ifndef CUBECLASS_HEADER
#define CUBECLASS_HEADER

#include <vector>
#include "utils.h"

#include "kdtree.h"
#include "index-list.h"

enum class CubeMethod {
  err = 0,
  cube = 1,
  lcube = 2
};

CubeMethod IntToCubeMethod(const int i);

class Cube {
protected:
  bool set_indirect = false;
  bool set_draw = false;

  void (Cube::*_Draw)(const int) = nullptr;

public:
  CubeMethod cubeMethod;

  int N;
  int p;
  double eps = 1e-12;

  IndexList* idx = nullptr;
  KDTree* tree = nullptr;
  std::vector<double> probabilities;

  std::vector<double> amat;
protected:
  std::vector<int> index;
  std::vector<int> neighbours;
  std::vector<double> distances;
  std::vector<double> bmat;
  std::vector<double> uvec;

public:
  std::vector<int> sample;

  // NO INIT CUBE
  Cube(
    const CubeMethod t_cubeMethod,
    const int t_N,
    const int t_p,
    const double t_eps
  );

  // DIRECT CUBE
  Cube(
    const double* t_probabilities,
    double* xxbalance,
    const int t_N,
    const int t_p,
    const double t_eps
  );
  // DIRECT LCUBE
  Cube(
    const double* t_probabilities,
    double* xxbalance,
    const int t_N,
    const int t_pbalance,
    const double t_eps,
    double* xxspread,
    const int t_pspread,
    const int bucketSize,
    const int method
  );

  ~Cube() {
    delete idx;
    delete tree;
  };

public:
  void AddUnitToSample(const int id);
  void EraseUnit(const int id);

  void InitIndirect();

protected:
  void Init(
    const double* t_probabilities,
    double* xxbalance,
    const int t_N,
    const int t_p,
    const double t_eps
  );

  int MaxSize();
  void RunUpdate();

  void Draw_cube(const int maxSize);
  void Draw_lcube(const int maxSize);

  void Draw(const int maxSize);
  void _RunFlight();
  void _RunLanding();

public:
  void RunFlight();
  void RunLanding();
  void Run();
};

#endif
