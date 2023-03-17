#ifndef CUBESTRATIFIEDCLASS_HEADER
#define CUBESTRATIFIEDCLASS_HEADER

#include <unordered_map>
#include <vector>
#include "CubeClass.h"

#include "kdtree.h"
#include "index-list.h"

class CubeStratified {
private:
  int N;
  int pb = 0;
  int ps = 0;
  double eps = 1e-12;

  CubeMethod cubeMethod;

  double* r_prob = nullptr; // NO DELETE
  double* r_xbalance = nullptr; // NO DELETE
  double* r_xspread = nullptr; // NO DELETE
  int* r_strata = nullptr; // NO DELETE

  int treeBucketSize = 40;
  int treeMethod = 2;
  KDTree* tree = nullptr;

  IndexList* idx = nullptr;

  std::unordered_map<int, int> stratumMap;

  Cube* cube = nullptr;

public:
  std::vector<int> sample;
private:
  std::vector<double> probabilities;
  std::vector<int> index;
  std::vector<int> stratumArr;
  double* spread = nullptr;

public:
  CubeStratified(
    double* t_prob,
    double* t_xbalance,
    int* t_strata,
    const int t_N,
    const int t_pb,
    const double t_eps
  );

  CubeStratified(
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
  );

  ~CubeStratified() {
    delete tree;
    delete idx;
    delete cube;
    delete[] spread;
  };

private:
  void Init(
    const int t_N,
    const int t_pb,
    const double t_eps,
    double* t_prob,
    double* t_xbalance,
    int* t_strata
  );

  void AddUnitToSample(const int id);
  void EraseUnit(const int id);

  void RunFlightPerStratum();
  void RunFlightOnFull();
  void RunLandingPerStratum();

public:
  void Run();

};

#endif
