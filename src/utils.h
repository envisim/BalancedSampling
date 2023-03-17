#ifndef BSUTILS_HEADER
#define BSUTILS_HEADER

inline bool Probability1(double p, double eps) {
  return (p) >= 1.0 - (eps);
};

inline bool Probability0(double p, double eps) {
  return (p) <= (eps);
};
inline bool ProbabilityInt(double p, double eps) {
  return (p) <= (eps) || (p) >= 1.0 - (eps);
};

inline int MatrixIdx(int row, int col, int pcols) {
  return (row) * (pcols) + (col);
}
inline int MatrixIdxRow(int row, int col, int pcols) {
  return (row) * (pcols) + (col);
}
inline int MatrixIdxCol(int row, int col, int nrows) {
  return (col) * (nrows) + (row);
}

#endif
