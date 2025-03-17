#include "svt_apply.h"

#include <cmath>
#include <functional>

#include "Rcpp.h"
#include "sparse_matrix.h"

using namespace Rcpp;

namespace smallcount {
namespace {

NumericVector getNzVals(const List svt_entry) {
    auto nz_vals = svt_entry[kSvtValInd];
    if (nz_vals == R_NilValue) {
        // Lacunar leaf (all ones)
        const auto nz_rows = as<IntegerVector>(svt_entry[kSvtRowInd]);
        return NumericVector(nz_rows.size(), 1);
    }
    return as<NumericVector>(nz_vals);
}

}  // namespace

List svtApply(Transformation transform, List old_svt, NumericVector mu) {
    const auto ncols = old_svt.size();
    List svt(ncols);
    size_t nz_index = 0;
    for (int i = 0; i < ncols; i++) {
        // NULL leaf (all zeros)
        if (old_svt[i] == R_NilValue) {
            continue;
        }

        const List old_entry = as<List>(old_svt[i]);
        const IntegerVector row_inds = as<IntegerVector>(old_entry[kSvtRowInd]);
        if (kSvtValInd != 0) {
            // Should never happen.
            stop("'nzvals' are stored at index %d of the SVT instead of 0",
                 kSvtValInd);
        }
        svt[i] = List::create(NumericVector(row_inds.size()), row_inds);

        // Apply transformation to all non-zero values.
        List entry = as<List>(svt[i]);
        const NumericVector old_nz_vals = getNzVals(old_entry);
        for (int j = 0; j < old_nz_vals.size(); j++) {
            NumericVector nz_vals = as<NumericVector>(entry[kSvtValInd]);
            nz_vals[j] = transform(old_nz_vals[j], mu[nz_index]);
            nz_index++;
        }
    }
    return svt;
}

}  // namespace smallcount
