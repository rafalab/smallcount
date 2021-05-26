#' use RcppArmadillo to parse 'sparse csv'
#' @import inline
#' @import Matrix
#' @param path character(1) path to uncompressed CSV with no header
#' @return list with component sp, a dgCMatrix instance
#' @examples
#' # superficial
#' if (!requireNamespace("bench")) stop("install bench package to run this example")
#' pa = system.file("extdata/tenx_subset-transpose.csv.gz", package="smallcount")
#' tf = tempfile()
#' file.copy(pa, tf)
#' file.rename(tf, ntf <- paste0(tf, ".gz"))
#' system(paste("gunzip", ntf)) # restores tf as gunzipped pa
#' tf2 = tempfile()
#' system(sprintf("sed -e '1,1d' %s | cut -d ',' -f 2- > %s", tf, tf2))
#' tim = bench::mark(mm <- parse_sparse_csv(tf2), iterations=2)
#' tim
#' dim(mm$sp)
#' mm$sp[11:16,1:6]
#' @export
parse_sparse_csv = function(path) {
 requireNamespace("Matrix") # for dgCMatrix
 g <- cxxfunction ( signature (vs="character"),
   plugin ="RcppArmadillo", body ='
   #include <armadillo>
   std::string v = Rcpp::as<std::string>(vs);
   arma::sp_mat D;
   D.load(v, arma::csv_ascii);
   return Rcpp::List::create(Rcpp::Named("sp")=D);
   ')
 g(path)
}
