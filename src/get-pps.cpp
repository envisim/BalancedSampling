#include <Rcpp.h>

// [[Rcpp::export(.getpps_cpp)]]
Rcpp::NumericVector getpps_cpp(
  Rcpp::NumericVector &x,
  int n
) {
  // Assume < 0 and equality has already been checked
  int N = x.length();
  Rcpp::NumericVector prob(N);
  int *index = new int[N];
  int one = 0;
  int oneplus = 0;
  int oneminus = 0;



  double xsum = 0.0;
  for (int i = 0; i < N; i++) {
    if (x[i] < 0.0)
      std::range_error("elements in x must be >= 0.0");

    xsum += x[i];
  }

  double temp = ((double) n) / xsum;
  for (int i = 0; i < N; i++) {
    prob[i] = x[i] * temp;

    if (prob[i] == 1.0) {
      one += 1;
    } else if (prob[i] > 1.0) {
      prob[i] = 1.0;
      oneplus += 1;
    } else {
      index[oneminus] = i;
      oneminus += 1;
    }
  }

  if (oneplus == 0) {
    delete[] index;
    return prob;
  }

  int onesum = one + oneplus;
  while (oneplus > 0) {
    xsum = 0.0;
    for (int i = 0; i < oneminus; i++)
      xsum += x[index[i]];

    temp = ((double) (n - onesum)) / xsum;

    int toneminus = oneminus;
    one = 0;
    oneplus = 0;
    oneminus = 0;


    for (int i = 0; i < toneminus; i++) {
      prob[index[i]] = x[index[i]] * temp;

      if (prob[index[i]] == 1.0) {
        one += 1;
      } else if (prob[index[i]] > 1.0) {
        prob[index[i]] = 1.0;
        oneplus += 1;
      } else {
        index[oneminus] = index[i];
        oneminus += 1;
      }
    }

    onesum += one + oneplus;
  }

  delete[] index;
  return prob;
}
