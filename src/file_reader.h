#ifndef SMALLCOUNT_FILE_READER_H_
#define SMALLCOUNT_FILE_READER_H_

#include <string>

#include "Rcpp.h"
#include "tenx_file_params.h"

using namespace Rcpp;

namespace smallcount {

// Static class to parse SparseMatrix objects.
class SparseMatrixFileReader {
   public:
    // Reads a SparseMatrix object from a file or directory.
    static SEXP read(const std::string &filepath, const std::string &rep,
                     const TenxFileParams &params);
};

}  // namespace smallcount

#endif
