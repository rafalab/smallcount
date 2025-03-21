// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>& Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppReadSparseMatrix
SEXP cppReadSparseMatrix(std::string sample, bool barcode_col_names,
                         bool id_row_names, std::string genome,
                         bool use_features_tsv);
RcppExport SEXP _smallcount_cppReadSparseMatrix(SEXP sampleSEXP,
                                                SEXP barcode_col_namesSEXP,
                                                SEXP id_row_namesSEXP,
                                                SEXP genomeSEXP,
                                                SEXP use_features_tsvSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter<std::string>::type sample(sampleSEXP);
    Rcpp::traits::input_parameter<bool>::type barcode_col_names(
        barcode_col_namesSEXP);
    Rcpp::traits::input_parameter<bool>::type id_row_names(id_row_namesSEXP);
    Rcpp::traits::input_parameter<std::string>::type genome(genomeSEXP);
    Rcpp::traits::input_parameter<bool>::type use_features_tsv(
        use_features_tsvSEXP);
    rcpp_result_gen = Rcpp::wrap(cppReadSparseMatrix(
        sample, barcode_col_names, id_row_names, genome, use_features_tsv));
    return rcpp_result_gen;
    END_RCPP
}
// cppPoissonDevianceTransformation
List cppPoissonDevianceTransformation(List svt, NumericVector mu);
RcppExport SEXP _smallcount_cppPoissonDevianceTransformation(SEXP svtSEXP,
                                                             SEXP muSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter<List>::type svt(svtSEXP);
    Rcpp::traits::input_parameter<NumericVector>::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(cppPoissonDevianceTransformation(svt, mu));
    return rcpp_result_gen;
    END_RCPP
}
// cppPoissonDispersionTransformation
List cppPoissonDispersionTransformation(List svt, NumericVector mu);
RcppExport SEXP _smallcount_cppPoissonDispersionTransformation(SEXP svtSEXP,
                                                               SEXP muSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter<List>::type svt(svtSEXP);
    Rcpp::traits::input_parameter<NumericVector>::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(cppPoissonDispersionTransformation(svt, mu));
    return rcpp_result_gen;
    END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_smallcount_cppReadSparseMatrix",
     (DL_FUNC)&_smallcount_cppReadSparseMatrix, 5},
    {"_smallcount_cppPoissonDevianceTransformation",
     (DL_FUNC)&_smallcount_cppPoissonDevianceTransformation, 2},
    {"_smallcount_cppPoissonDispersionTransformation",
     (DL_FUNC)&_smallcount_cppPoissonDispersionTransformation, 2},
    {NULL, NULL, 0}};

RcppExport void R_init_smallcount(DllInfo* dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
