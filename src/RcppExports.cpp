// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cube_stratified_cpp
Rcpp::IntegerVector cube_stratified_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, Rcpp::IntegerVector& strata, double eps);
RcppExport SEXP _BalancedSampling_cube_stratified_cpp(SEXP probSEXP, SEXP xSEXP, SEXP strataSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cube_stratified_cpp(prob, x, strata, eps));
    return rcpp_result_gen;
END_RCPP
}
// cube_cpp
Rcpp::IntegerVector cube_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, double eps);
RcppExport SEXP _BalancedSampling_cube_cpp(SEXP probSEXP, SEXP xSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cube_cpp(prob, x, eps));
    return rcpp_result_gen;
END_RCPP
}
// cube_fast_cpp
Rcpp::IntegerVector cube_fast_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, double eps);
RcppExport SEXP _BalancedSampling_cube_fast_cpp(SEXP probSEXP, SEXP xSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cube_fast_cpp(prob, x, eps));
    return rcpp_result_gen;
END_RCPP
}
// hlpm2_cpp
Rcpp::IntegerMatrix hlpm2_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, Rcpp::IntegerVector& sizes, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_hlpm2_cpp(SEXP probSEXP, SEXP xSEXP, SEXP sizesSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(hlpm2_cpp(prob, x, sizes, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// lcps
Rcpp::NumericVector lcps(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x);
RcppExport SEXP _BalancedSampling_lcps(SEXP probSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(lcps(prob, x));
    return rcpp_result_gen;
END_RCPP
}
// lcube_stratified_cpp
Rcpp::IntegerVector lcube_stratified_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& xbal, Rcpp::NumericMatrix& xspread, Rcpp::IntegerVector& strata, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_lcube_stratified_cpp(SEXP probSEXP, SEXP xbalSEXP, SEXP xspreadSEXP, SEXP strataSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type xbal(xbalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type xspread(xspreadSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(lcube_stratified_cpp(prob, xbal, xspread, strata, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// lcube_cpp
Rcpp::IntegerVector lcube_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& xbal, Rcpp::NumericMatrix& xspread, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_lcube_cpp(SEXP probSEXP, SEXP xbalSEXP, SEXP xspreadSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type xbal(xbalSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type xspread(xspreadSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(lcube_cpp(prob, xbal, xspread, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// lpm1_search_cpp
Rcpp::IntegerVector lpm1_search_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_lpm1_search_cpp(SEXP probSEXP, SEXP xSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(lpm1_search_cpp(prob, x, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// lpm1_cpp
Rcpp::IntegerVector lpm1_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_lpm1_cpp(SEXP probSEXP, SEXP xSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(lpm1_cpp(prob, x, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// lpm2_int_cpp
Rcpp::IntegerVector lpm2_int_cpp(int n, Rcpp::NumericMatrix& x, int bucketSize, int method);
RcppExport SEXP _BalancedSampling_lpm2_int_cpp(SEXP nSEXP, SEXP xSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(lpm2_int_cpp(n, x, bucketSize, method));
    return rcpp_result_gen;
END_RCPP
}
// lpm2_cpp
Rcpp::IntegerVector lpm2_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_lpm2_cpp(SEXP probSEXP, SEXP xSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(lpm2_cpp(prob, x, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// rpm_cpp
Rcpp::IntegerVector rpm_cpp(Rcpp::NumericVector& prob, double eps);
RcppExport SEXP _BalancedSampling_rpm_cpp(SEXP probSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(rpm_cpp(prob, eps));
    return rcpp_result_gen;
END_RCPP
}
// sb_sumofsquares_cpp
double sb_sumofsquares_cpp(Rcpp::NumericMatrix& x, Rcpp::IntegerVector& sample, int bucketSize, int method);
RcppExport SEXP _BalancedSampling_sb_sumofsquares_cpp(SEXP xSEXP, SEXP sampleSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sb_sumofsquares_cpp(x, sample, bucketSize, method));
    return rcpp_result_gen;
END_RCPP
}
// sb_voronoi_cpp
double sb_voronoi_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, Rcpp::IntegerVector& sample, int bucketSize, int method);
RcppExport SEXP _BalancedSampling_sb_voronoi_cpp(SEXP probSEXP, SEXP xSEXP, SEXP sampleSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(sb_voronoi_cpp(prob, x, sample, bucketSize, method));
    return rcpp_result_gen;
END_RCPP
}
// scps_coord_cpp
Rcpp::IntegerVector scps_coord_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, Rcpp::NumericVector& random, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_scps_coord_cpp(SEXP probSEXP, SEXP xSEXP, SEXP randomSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type random(randomSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(scps_coord_cpp(prob, x, random, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// scps_cpp
Rcpp::IntegerVector scps_cpp(Rcpp::NumericVector& prob, Rcpp::NumericMatrix& x, int bucketSize, int method, double eps);
RcppExport SEXP _BalancedSampling_scps_cpp(SEXP probSEXP, SEXP xSEXP, SEXP bucketSizeSEXP, SEXP methodSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(scps_cpp(prob, x, bucketSize, method, eps));
    return rcpp_result_gen;
END_RCPP
}
// spm_cpp
Rcpp::IntegerVector spm_cpp(Rcpp::NumericVector& prob, double eps);
RcppExport SEXP _BalancedSampling_spm_cpp(SEXP probSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(spm_cpp(prob, eps));
    return rcpp_result_gen;
END_RCPP
}
// vsb
double vsb(NumericVector probs, NumericVector ys, NumericMatrix xs);
RcppExport SEXP _BalancedSampling_vsb(SEXP probsSEXP, SEXP ysSEXP, SEXP xsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xs(xsSEXP);
    rcpp_result_gen = Rcpp::wrap(vsb(probs, ys, xs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BalancedSampling_cube_stratified_cpp", (DL_FUNC) &_BalancedSampling_cube_stratified_cpp, 4},
    {"_BalancedSampling_cube_cpp", (DL_FUNC) &_BalancedSampling_cube_cpp, 3},
    {"_BalancedSampling_cube_fast_cpp", (DL_FUNC) &_BalancedSampling_cube_fast_cpp, 3},
    {"_BalancedSampling_hlpm2_cpp", (DL_FUNC) &_BalancedSampling_hlpm2_cpp, 6},
    {"_BalancedSampling_lcps", (DL_FUNC) &_BalancedSampling_lcps, 2},
    {"_BalancedSampling_lcube_stratified_cpp", (DL_FUNC) &_BalancedSampling_lcube_stratified_cpp, 7},
    {"_BalancedSampling_lcube_cpp", (DL_FUNC) &_BalancedSampling_lcube_cpp, 6},
    {"_BalancedSampling_lpm1_search_cpp", (DL_FUNC) &_BalancedSampling_lpm1_search_cpp, 5},
    {"_BalancedSampling_lpm1_cpp", (DL_FUNC) &_BalancedSampling_lpm1_cpp, 5},
    {"_BalancedSampling_lpm2_int_cpp", (DL_FUNC) &_BalancedSampling_lpm2_int_cpp, 4},
    {"_BalancedSampling_lpm2_cpp", (DL_FUNC) &_BalancedSampling_lpm2_cpp, 5},
    {"_BalancedSampling_rpm_cpp", (DL_FUNC) &_BalancedSampling_rpm_cpp, 2},
    {"_BalancedSampling_sb_sumofsquares_cpp", (DL_FUNC) &_BalancedSampling_sb_sumofsquares_cpp, 4},
    {"_BalancedSampling_sb_voronoi_cpp", (DL_FUNC) &_BalancedSampling_sb_voronoi_cpp, 5},
    {"_BalancedSampling_scps_coord_cpp", (DL_FUNC) &_BalancedSampling_scps_coord_cpp, 6},
    {"_BalancedSampling_scps_cpp", (DL_FUNC) &_BalancedSampling_scps_cpp, 5},
    {"_BalancedSampling_spm_cpp", (DL_FUNC) &_BalancedSampling_spm_cpp, 2},
    {"_BalancedSampling_vsb", (DL_FUNC) &_BalancedSampling_vsb, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_BalancedSampling(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
