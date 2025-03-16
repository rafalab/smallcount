#include "sparse_matrix.h"

#include <memory>
#include <string>
#include <vector>

#include "Rcpp.h"

using namespace Rcpp;

namespace smallcount {
namespace {

static constexpr char kSvtSparseMatrix[] = "SVT_SparseMatrix";

static constexpr char kSvt[] = "SVT";
static constexpr char kSvtVersion[] = ".svt_version";
static constexpr char kType[] = "type";
static constexpr char kDim[] = "dim";
static constexpr char kDimNames[] = "dimnames";

static constexpr int kVersionNum = 1;
static constexpr char kInteger[] = "integer";

List createDimNamesList(std::vector<std::string> row_names,
                        std::vector<std::string> col_names) {
    List dim_names = List(2);
    if (!row_names.empty()) {
        dim_names[0] = wrap(std::move(row_names));
    }
    if (!col_names.empty()) {
        dim_names[1] = wrap(std::move(col_names));
    }
    return dim_names;
}

void sortRowIndices(SvtEntry &col) {
    // Check if already sorted.
    if (std::is_sorted(col[kSvtRowInd].begin(), col[kSvtRowInd].end())) {
        return;
    }

    // Check if reverse sorted.
    if (std::is_sorted(col[kSvtRowInd].begin(), col[kSvtRowInd].end(),
                       std::greater<int>())) {
        std::reverse(col[kSvtRowInd].begin(), col[kSvtRowInd].end());
        std::reverse(col[kSvtValInd].begin(), col[kSvtValInd].end());
        return;
    }

    // Create vector of indices sorted by row.
    std::vector<size_t> indices(col[kSvtRowInd].size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return col[kSvtRowInd][i] < col[kSvtRowInd][j];
    });

    // Sort parallel arrays by row index.
    std::vector<int> sorted_rows(col[kSvtRowInd].size());
    std::vector<int> sorted_vals(col[kSvtValInd].size());
    for (size_t i = 0; i < indices.size(); ++i) {
        sorted_rows[i] = col[kSvtRowInd][indices[i]];
        sorted_vals[i] = col[kSvtValInd][indices[i]];
    }
    col[kSvtRowInd] = std::move(sorted_rows);
    col[kSvtValInd] = std::move(sorted_vals);
}

}  // namespace

List SvtSparseMatrix::createSvtList(Svt svt) {
    for (auto &col : svt) {
        sortRowIndices(col);
    }
    const int ncol = svt.size();
    List svt_list = List(ncol);
    for (int i = 0; i < ncol; i++) {
        if (!svt[i][0].empty()) {
            svt_list[i] = wrap(std::move(svt[i]));
        }
    }
    return svt_list;
}

SEXP SvtSparseMatrix::toRcpp() {
    S4 obj(kSvtSparseMatrix);
    obj.slot(kSvt) = R_NilValue;
    if (metadata.nval != 0) {
        obj.slot(kSvt) = createSvtList(std::move(svt));
    }
    obj.slot(kDim) = IntegerVector({metadata.nrow, metadata.ncol});
    obj.slot(kDimNames) = createDimNamesList(std::move(metadata.row_names),
                                             std::move(metadata.col_names));
    obj.slot(kType) = kInteger;
    obj.slot(kSvtVersion) = kVersionNum;
    return obj;
}

}  // namespace smallcount
