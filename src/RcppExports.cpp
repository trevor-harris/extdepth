// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// depth
NumericVector depth(NumericVector f, NumericMatrix fmat);
RcppExport SEXP _extdepth_depth(SEXP fSEXP, SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(depth(f, fmat));
    return rcpp_result_gen;
END_RCPP
}
// dCDF
NumericVector dCDF(NumericVector depths);
RcppExport SEXP _extdepth_dCDF(SEXP depthsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type depths(depthsSEXP);
    rcpp_result_gen = Rcpp::wrap(dCDF(depths));
    return rcpp_result_gen;
END_RCPP
}
// ed_compare
int ed_compare(NumericVector f1_cdf, NumericVector f2_cdf);
RcppExport SEXP _extdepth_ed_compare(SEXP f1_cdfSEXP, SEXP f2_cdfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f1_cdf(f1_cdfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f2_cdf(f2_cdfSEXP);
    rcpp_result_gen = Rcpp::wrap(ed_compare(f1_cdf, f2_cdf));
    return rcpp_result_gen;
END_RCPP
}
// ED
NumericVector ED(NumericMatrix fmat);
RcppExport SEXP _extdepth_ED(SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(ED(fmat));
    return rcpp_result_gen;
END_RCPP
}
// fast_depth
NumericVector fast_depth(NumericVector f, NumericMatrix fmat);
RcppExport SEXP _extdepth_fast_depth(SEXP fSEXP, SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_depth(f, fmat));
    return rcpp_result_gen;
END_RCPP
}
// dCDF_r
double dCDF_r(NumericVector depths, double r);
RcppExport SEXP _extdepth_dCDF_r(SEXP depthsSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type depths(depthsSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(dCDF_r(depths, r));
    return rcpp_result_gen;
END_RCPP
}
// fast_compare
int fast_compare(NumericVector d1, NumericVector d2);
RcppExport SEXP _extdepth_fast_compare(SEXP d1SEXP, SEXP d2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d2(d2SEXP);
    rcpp_result_gen = Rcpp::wrap(fast_compare(d1, d2));
    return rcpp_result_gen;
END_RCPP
}
// fast_ED
NumericVector fast_ED(NumericMatrix fmat);
RcppExport SEXP _extdepth_fast_ED(SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_ED(fmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_extdepth_depth", (DL_FUNC) &_extdepth_depth, 2},
    {"_extdepth_dCDF", (DL_FUNC) &_extdepth_dCDF, 1},
    {"_extdepth_ed_compare", (DL_FUNC) &_extdepth_ed_compare, 2},
    {"_extdepth_ED", (DL_FUNC) &_extdepth_ED, 1},
    {"_extdepth_fast_depth", (DL_FUNC) &_extdepth_fast_depth, 2},
    {"_extdepth_dCDF_r", (DL_FUNC) &_extdepth_dCDF_r, 2},
    {"_extdepth_fast_compare", (DL_FUNC) &_extdepth_fast_compare, 2},
    {"_extdepth_fast_ED", (DL_FUNC) &_extdepth_fast_ED, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_extdepth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
