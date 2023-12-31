// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// scale
NumericMatrix scale(const NumericMatrix& X);
RcppExport SEXP _dummyscoring_scale(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(scale(X));
    return rcpp_result_gen;
END_RCPP
}
// sparsecor
NumericMatrix sparsecor(NumericMatrix& X, NumericMatrix& Y);
RcppExport SEXP _dummyscoring_sparsecor(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(sparsecor(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// score_matrix
NumericMatrix score_matrix(arma::mat& X, arma::mat& P, arma::mat& Uinv, arma::mat& R, String method, int nthreads, bool showprogress);
RcppExport SEXP _dummyscoring_score_matrix(SEXP XSEXP, SEXP PSEXP, SEXP UinvSEXP, SEXP RSEXP, SEXP methodSEXP, SEXP nthreadsSEXP, SEXP showprogressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Uinv(UinvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type showprogress(showprogressSEXP);
    rcpp_result_gen = Rcpp::wrap(score_matrix(X, P, Uinv, R, method, nthreads, showprogress));
    return rcpp_result_gen;
END_RCPP
}
// dummy_extension
List dummy_extension(NumericMatrix& X, NumericMatrix& X_v, NumericMatrix& P_vf, arma::mat& R_vv);
RcppExport SEXP _dummyscoring_dummy_extension(SEXP XSEXP, SEXP X_vSEXP, SEXP P_vfSEXP, SEXP R_vvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type X_v(X_vSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type P_vf(P_vfSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type R_vv(R_vvSEXP);
    rcpp_result_gen = Rcpp::wrap(dummy_extension(X, X_v, P_vf, R_vv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dummyscoring_scale", (DL_FUNC) &_dummyscoring_scale, 1},
    {"_dummyscoring_sparsecor", (DL_FUNC) &_dummyscoring_sparsecor, 2},
    {"_dummyscoring_score_matrix", (DL_FUNC) &_dummyscoring_score_matrix, 7},
    {"_dummyscoring_dummy_extension", (DL_FUNC) &_dummyscoring_dummy_extension, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dummyscoring(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
