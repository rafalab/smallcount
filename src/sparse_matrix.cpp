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
static constexpr char kType[] = "type";
static constexpr char kInteger[] = "integer";

static constexpr char kDim[] = "dim";
static constexpr char kDimNames[] = "dimnames";

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

}  // namespace

List SvtSparseMatrix::createSvtList(Svt svt) {
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
    return obj;
}

}  // namespace smallcount
