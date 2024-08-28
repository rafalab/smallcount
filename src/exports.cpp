#include <string>

#include "Rcpp.h"
#include "file_reader.h"

using namespace Rcpp;

// Reads a SparseMatrix R object from a file using the internal representation
// specified (either "coo" or "svt").
// [[Rcpp::export]]
SEXP cppReadSparseMatrix(std::string filepath, std::string rep) {
    return smallcount::SparseMatrixFileReader::read(filepath, rep);
}
