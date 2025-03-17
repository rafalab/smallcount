#ifndef SMALLCOUNT_SVT_APPLY_H_
#define SMALLCOUNT_SVT_APPLY_H_

#include "Rcpp.h"

using namespace Rcpp;

namespace smallcount {

using Transformation = std::function<double(double, double)>;

// Performs `nzvals <- transform(nzvals, mu)`
List svtApply(Transformation transform, List old_svt, NumericVector mu);

}  // namespace smallcount

#endif  // SMALLCOUNT_SVT_APPLY_H_
