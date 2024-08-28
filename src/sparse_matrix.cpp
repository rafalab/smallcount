#include "sparse_matrix.h"

#include <memory>
#include <string>

#include "Rcpp.h"

using namespace Rcpp;

namespace smallcount {
namespace {

static constexpr char kCooSparseMatrix[] = "COO_SparseMatrix";
static constexpr char kNonZeroCoords[] = "nzcoo";
static constexpr char kNonZeroData[] = "nzdata";

static constexpr char kSvtSparseMatrix[] = "SVT_SparseMatrix";
static constexpr char kSvt[] = "SVT";
static constexpr char kType[] = "type";
static constexpr char kInteger[] = "integer";
static constexpr int kSvtRowInd = 0;  // Index of row information in SVT entry
static constexpr int kSvtValInd = 1;  // Index of value information in SVT entry

static constexpr char kDim[] = "dim";

}  // namespace

std::unique_ptr<SparseMatrix> SparseMatrix::create(const std::string &rep) {
    if (rep == kCooRep) {
        return std::make_unique<CooSparseMatrix>();
    }
    if (rep == kSvtRep) {
        return std::make_unique<SvtSparseMatrix>();
    }
    return nullptr;
}

// -----------------------------------------------------------------------------
// COO Sparse Matrix
// -----------------------------------------------------------------------------
void CooSparseMatrix::init(const MatrixMetadata &m) {
    metadata = m;
    rows.reserve(m.nval);
    cols.reserve(m.nval);
    vals.reserve(m.nval);
}

void CooSparseMatrix::addEntry(const MatrixData &entry) {
    rows.push_back(entry.row);
    cols.push_back(entry.col);
    vals.push_back(entry.val);
}

SEXP CooSparseMatrix::toRcpp() {
    S4 obj(kCooSparseMatrix);
    obj.slot(kNonZeroCoords) =
        createCoordsMatrix(std::move(rows), std::move(cols));
    obj.slot(kNonZeroData) = wrap(std::move(vals));
    obj.slot(kDim) = IntegerVector({nrow(), ncol()});
    return obj;
}

IntegerMatrix CooSparseMatrix::createCoordsMatrix(std::vector<int> rows,
                                                  std::vector<int> cols) {
    IntegerMatrix coords(rows.size(), 2);
    std::move(rows.begin(), rows.end(), coords.begin());
    std::move(cols.begin(), cols.end(), coords.begin() + rows.size());
    return coords;
}

// -----------------------------------------------------------------------------
// SVT Sparse Matrix
// -----------------------------------------------------------------------------
void SvtSparseMatrix::init(const MatrixMetadata &m) {
    metadata = m;
    svt = Svt(m.ncol);
}

void SvtSparseMatrix::addEntry(const MatrixData &entry) {
    // Shift row and column by -1 for 0-based indexing.
    const int col = entry.col - 1;
    if (svt[col].empty()) {
        svt[col] = SvtEntry(2);
    }
    svt[col][kSvtRowInd].push_back(entry.row - 1);
    svt[col][kSvtValInd].push_back(entry.val);
}

SEXP SvtSparseMatrix::toRcpp() {
    S4 obj(kSvtSparseMatrix);
    obj.slot(kSvt) = R_NilValue;
    if (nval() != 0) {
        obj.slot(kSvt) = createSvtList(std::move(svt));
    }
    obj.slot(kDim) = IntegerVector({nrow(), ncol()});
    obj.slot(kType) = kInteger;
    return obj;
}

List SvtSparseMatrix::createSvtList(Svt svt) {
    const int ncol = svt.size();
    List svt_list = List(ncol);
    for (int i = 0; i < ncol; i++) {
        if (!svt[i].empty()) {
            svt_list[i] = wrap(std::move(svt[i]));
        }
    }
    return svt_list;
}

}  // namespace smallcount
