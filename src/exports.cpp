#include <string>

#include "Rcpp.h"
#include "file_reader.h"
#include "tenx_file_params.h"

using namespace Rcpp;

// Reads a SparseMatrix object from a file or directory.
// [[Rcpp::export]]
SEXP cppReadSparseMatrix(std::string sample, bool barcode_col_names,
                         bool id_row_names, std::string genome,
                         bool use_features_tsv, std::string rep) {
    smallcount::TenxFileParams file_params;
    file_params.use_barcode_col_names = barcode_col_names;
    file_params.use_id_row_names = id_row_names;
    if (!genome.empty()) {
        file_params.genome.emplace(genome);
    }
    file_params.use_features_tsv = use_features_tsv;
    return smallcount::SparseMatrixFileReader::read(sample, rep, file_params);
}
