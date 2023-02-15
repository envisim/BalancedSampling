#include "uniform.h"

//**********************************************
// Thx to R
//**********************************************

double stduniform() {
  double u;
  do {u = R::unif_rand();} while (u < 0.0 || u >= 1.0);
  return u;
}
