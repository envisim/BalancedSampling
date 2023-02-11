#include <Rcpp.h>

//**********************************************
// Author: Wilmer Prentius
// Last edit: 2023-02-08
// Licence: GPL (>=2)
//**********************************************

// struct Object {
//   int probability;
//   unsigned int index;
//   Object() : probability(0), index(0) {}
// };

double euclideanDistance(double *x, int id1, int id2, int cols) {
  double dist = 0.0;
  double *u1 = x + id1 * cols;
  double *u2 = x + id2 * cols;

  for (int k = 0; k < cols; k++) {
    double temp = *(u1 + k) - *(u2 + k);
    dist += temp * temp;
  }

  return dist;
}

int randn(double u, int n) {
  return (int)((double)n * u);
}

// [[Rcpp::export]]
Rcpp::IntegerVector lpm2int(int n, Rcpp::NumericMatrix &x) {
  int N = x.ncol();
  int P = x.nrow();
  int unresolvedObjects = N;
  double *xx = REAL(x);

  int *probability = new int[N];
  int *idx = new int[N];
  int *neighbours = new int[N];

  for (int i = 0; i < N; i++) {
    probability[i] = n;
    idx[i] = i;
  }

  Rcpp::NumericVector rand1 = Rcpp::runif(N, 0.0, 1.0);
  Rcpp::NumericVector rand2 = Rcpp::runif(N, 0.0, 1.0);
  // Rcpp::NumericVector rand3 = R::runif(N);

  while (unresolvedObjects > 1) {
    int u1 = randn(rand2[unresolvedObjects-1], unresolvedObjects);
    // int u1 = randn(R::runif(0.0, 1.0), unresolvedObjects);
    int idx1 = idx[u1];

    double mindist = DBL_MAX;
    int len = 0;

    for (int i = 0; i < unresolvedObjects; i++) {
      if (u1 == i)
        continue;

      double dist = euclideanDistance(xx, idx1, idx[i], P);

      if (dist < mindist) {
        mindist = dist;
        neighbours[0] = i;
        len = 1;
      } else if (dist == mindist) {
        neighbours[len] = i;
        len += 1;
      }
    }

    // int u2 = len == 1 ? neighbours[0] : randn(rand3[unresolvedObjects], len);
    int u2 = len == 1 ? neighbours[0] : randn(R::runif(0.0, 1.0), len);
    int idx2 = idx[u2];

    int p1 = probability[idx1];
    int p2 = probability[idx2];
    int psum = p1 + p2;
    double u = rand1[unresolvedObjects-1];
    // double u = R::runif(0.0, 1.0);

    if (psum > N) {
      if (N - p2 > randn(u, (N << 1) - psum)) {
        probability[idx1] = N;
        probability[idx2] = psum - N;
      } else {
        probability[idx1] = psum - N;
        probability[idx2] = N;
      }
    } else {
      if (p2 > randn(u, psum)) {
        probability[idx1] = 0;
        probability[idx2] = psum;
      } else {
        probability[idx1] = psum;
        probability[idx2] = 0;
      }
    }

    if (probability[idx1] == 0 || probability[idx1] == N) {
      unresolvedObjects -= 1;
      int temp = idx[unresolvedObjects];
      idx[unresolvedObjects] = idx1;
      idx[u1] = temp;

      if (unresolvedObjects == u2)
        u2 = u1;
    }

    if (probability[idx2] == 0 || probability[idx2] == N) {
      unresolvedObjects -= 1;
      int temp = idx[unresolvedObjects];
      idx[unresolvedObjects] = idx2;
      idx[u2] = temp;
    }
  }

  Rcpp::IntegerVector s(n);

  for (int i = 0, j = 0; i < N && j < n; i++) {
    if (probability[i] == N) {
      s[j] = i + 1;
      j += 1;
    }
  }

  delete[] probability;
  delete[] idx;
  delete[] neighbours;

  return s;
}

