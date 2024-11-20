#include <string>

#include "Rcpp.h"
#include "file_params.h"
#include "file_reader.h"

using namespace Rcpp;

// Reads a SparseMatrix object from a file or directory.
// [[Rcpp::export]]
SEXP cppReadSparseMatrix(std::string sample, bool col_names,
                         std::string row_names, std::string genome,
                         std::string rep) {
    smallcount::FileParams file_params;
    file_params.use_barcode_col_names = col_names;
    file_params.use_id_row_names = (row_names == "id");
    if (!genome.empty()) {
        file_params.genome.emplace(genome);
    }
    return smallcount::SparseMatrixFileReader::read(sample, rep, file_params);
}
