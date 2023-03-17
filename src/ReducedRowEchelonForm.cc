#include "ReducedRowEchelonForm.h"
#include "utils.h"

void ReducedRowEchelonForm(
  double *mat,
  const int rowCount,
  const int colCount
) {
  int lead = 0;

  for (int r = 0; r < rowCount; r++) {
    if (colCount <= lead)
      return;

    int i = r;

    while (mat[MatrixIdx(i, lead, colCount)] == 0) {
      i += 1;

      if (i == rowCount) {
        i = r;
        lead += 1;

        if (colCount == lead)
          return;
      }
    }

    double *br = mat + MatrixIdx(r, 0, colCount);

    if (i != r) {
      double *bi = mat + MatrixIdx(i, 0, colCount);

      for (int k = 0; k < colCount; k++) {
        double temp = bi[k];
        bi[k] = br[k];
        br[k] = temp;
      }
    }

    if (br[lead] != 0.0) {
      double temp = br[lead];
      for (int k = 0; k < colCount; k++) {
        br[k] /= temp;
      }
    }

    for (int j = 0; j < rowCount; j++) {
      if (j == r)
        continue;

      double temp = mat[MatrixIdx(j, lead, colCount)];
      for (int k = 0; k < colCount; k++) {
        mat[MatrixIdx(j, k, colCount)] -= br[k] * temp;
      }
    }

    lead += 1;
  }

  return;
}
