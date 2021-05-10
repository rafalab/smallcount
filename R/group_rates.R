#' Rowwise rates for groups
#'
#' @param y A tgCMatrix sparse Matrix.
#' @param g A factor defining the group for each column.
#'
#' @export
#'
group_rates <- function(y, g){

  if(!is(y, "dgCMatrix")) stop("y must be class dgCMatrix")

  if(!is.factor(g)){
    warning("Coercing g into a factor")
    g <- as.factor(g)
  }

  js <- as.numeric(g)

  rowsums <- matrix(0, nrow(y),  length(n))
  colsums <- vector("numeric", length(n))

  for(j in 1:ncol(y)){
    ind <- (y@p[j]+1):y@p[j+1]
    real_ind <- y@i[ind] + 1
    k <- js[j]
    x <- y@x[ind]
    rowsums[real_ind, k] <- rowsums[real_ind, k] + x
    colsums[k] <- colsums[k] + sum(x)
  }

  rowsums <- sweep(rowsums, 2, colsums, FUN = "/")
  colnames(rowsums) <- levels(g)
  rownames(rowsums) <- rownames(y)

  return(rowsums)

}

