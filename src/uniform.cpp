#include "uniform.h"

//**********************************************
// Thx to R
//**********************************************

double stduniform() {
  double u;
  do {u = R::unif_rand();} while (u < 0.0 || u >= 1.0);
  return u;
}

int intuniform(int N) {
  if (N == 1 || N == 0)
    return 0;

  return (int)((double)N * stduniform());
}

int intuniform(double N) {
  return (int)(N * stduniform());
}
