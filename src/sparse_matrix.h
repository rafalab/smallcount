#ifndef SMALLCOUNT_SPARSE_MATRIX_H_
#define SMALLCOUNT_SPARSE_MATRIX_H_

#include <memory>
#include <string>
#include <vector>

#include "Rcpp.h"

using namespace Rcpp;

namespace smallcount {

static constexpr char kCooRep[] = "coo";
static constexpr char kSvtRep[] = "svt";

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

// Non-zero entry in a sparse matrix.
struct MatrixData {
    int row;  // Row index (1-based index, for consistency with R)
    int col;  // Column index (1-based index, for consistency with R)
    int val;  // Value
};

// Abstract sparse matrix struct.
struct SparseMatrix {
    // Initializes the matrix with the given metadata.
    virtual void init(const MatrixMetadata &metadata) = 0;
    // Adds a non-zero data entry to the matrix.
    virtual void addEntry(const MatrixData &entry) = 0;
    // Consumes the matrix and converts it to an S4 object.
    virtual SEXP toRcpp() = 0;

    // Returns a sparse matrix with the specified internal representation.
    // Only "coo" and "svt" are currently supported.
    static std::unique_ptr<SparseMatrix> create(const std::string &rep);

    int nrow() const { return metadata.nrow; }
    int ncol() const { return metadata.ncol; }
    size_t nval() const { return metadata.nval; }

    virtual ~SparseMatrix() = default;

   protected:
    // Metadata for the number of rows, columns, and non-zero values.
    // Default initialization is an empty 1 x 1 matrix.
    MatrixMetadata metadata = {.nrow = 1, .ncol = 1, .nval = 0};

    // Checks that the matrix entry is not out of bounds.
    void checkValid(const MatrixData &entry);
};

// COO representation of a sparse matrix.
struct CooSparseMatrix : public SparseMatrix {
    // Row indices of non-zero entries.
    std::vector<int> rows{};
    // Column indices of non-zero entries.
    std::vector<int> cols{};
    // Non-zero values.
    std::vector<int> vals{};

    void init(const MatrixMetadata &metadata) override;
    void addEntry(const MatrixData &entry) override;
    SEXP toRcpp() override;

   private:
    // Merges the C++ coordinate vectors into an Rcpp IntegerMatrix.
    static IntegerMatrix createCoordsMatrix(std::vector<int> rows,
                                            std::vector<int> cols);
};

// SVT representation of a sparse matrix.
struct SvtSparseMatrix : public SparseMatrix {
    // Sparse vector tree.
    Svt svt{};

    void init(const MatrixMetadata &metadata) override;
    void addEntry(const MatrixData &entry) override;
    SEXP toRcpp() override;

   private:
    // Converts the C++ SVT matrix to an Rcpp List.
    static List createSvtList(Svt svt);
};

}  // namespace smallcount

#endif
