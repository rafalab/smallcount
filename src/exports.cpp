#include <string>

#include "Rcpp.h"
#include "file_reader.h"
#include "svt_apply.h"
#include "tenx_file_params.h"

using namespace Rcpp;
using smallcount::Transformation;

// Reads a SparseMatrix object from a file or directory.
// [[Rcpp::export]]
SEXP cppReadSparseMatrix(std::string sample, bool barcode_col_names,
                         bool id_row_names, std::string genome,
                         bool use_features_tsv) {
    smallcount::TenxFileParams file_params;
    file_params.use_barcode_col_names = barcode_col_names;
    file_params.use_id_row_names = id_row_names;
    if (!genome.empty()) {
        file_params.genome.emplace(genome);
    }
    file_params.use_features_tsv = use_features_tsv;
    return smallcount::SparseMatrixFileReader::read(sample, file_params);
}

// Performs `nzvals <- nzvals * log(nzvals / mu)`
// [[Rcpp::export]]
List cppPoissonDevianceTransformation(List svt, NumericVector mu) {
    Transformation dev = [](double y, double mu) { return y * log(y / mu); };
    return smallcount::svtApply(dev, svt, mu);
}

// Performs `nzvals <- nzvals^2 / mu`
// [[Rcpp::export]]
List cppPoissonDispersionTransformation(List svt, NumericVector mu) {
    Transformation disp = [](double y, double mu) { return y * y / mu; };
    return smallcount::svtApply(disp, svt, mu);
}
