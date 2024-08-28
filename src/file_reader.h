#ifndef SMALLCOUNT_FILE_READER_H_
#define SMALLCOUNT_FILE_READER_H_

#include <string>

#include "Rcpp.h"

using namespace Rcpp;

namespace smallcount {

// Static class to parse a sparse matrix from a file.
class SparseMatrixFileReader {
   public:
    // Reads a sparse matrix from a file using the internal representation
    // specified (either "coo" or "svt").
    static SEXP read(const std::string &filepath, const std::string &rep);
};

}  // namespace smallcount

#endif
