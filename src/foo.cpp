// include RcppArmadillo and Rcpp
#include "RcppArmadillo.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List parse_sparse_csv_impl(SEXP fname) {
using namespace Rcpp;
std::string v = Rcpp::as<std::string>(fname);
arma::sp_mat D;
D.load(v, arma::csv_ascii);
return Rcpp::List::create(Rcpp::Named("sp")=D);
}
