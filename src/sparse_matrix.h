#ifndef SMALLCOUNT_SPARSE_MATRIX_H_
#define SMALLCOUNT_SPARSE_MATRIX_H_

#include <memory>
#include <string>
#include <vector>

#include "Rcpp.h"

using namespace Rcpp;

namespace smallcount {

static constexpr int kSvtValInd = 0;  // Index of value information in SVT entry
static constexpr int kSvtRowInd = 1;  // Index of row information in SVT entry

// Pair of parallel arrays containing the row and value, respectively, of each
// non-zero entry in a given column.
using SvtEntry = std::vector<std::vector<int>>;
// Vector of SVT entries for each column.
using Svt = std::vector<SvtEntry>;

// Information about a sparse matrix.
struct MatrixMetadata {
    int nrow;     // Number of rows
    int ncol;     // Number of columns
    size_t nval;  // Number of non-zero values

    std::vector<std::string> row_names;  // Row names
    std::vector<std::string> col_names;  // Column names
};

// SVT representation of a sparse matrix.
struct SvtSparseMatrix {
    // Construct from pre-built Svt and metadata
    SvtSparseMatrix(Svt svt, MatrixMetadata metadata)
        : svt(std::move(svt)), metadata(std::move(metadata)) {}

    // Consumes the matrix and converts it to an S4 object.
    SEXP toRcpp();

   private:
    // Sparse vector tree.
    Svt svt{};
    // Metadata for the number of rows, columns, and non-zero values.
    MatrixMetadata metadata;

    // Converts the C++ SVT matrix to an Rcpp List.
    static List createSvtList(Svt svt);
};

}  // namespace smallcount

#endif
