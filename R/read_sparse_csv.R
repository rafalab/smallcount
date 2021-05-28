#' Read a sparse delimited file in a dgCMatrix
#' @note The functions are prototypes. They need to be ported to C++
#' The functions reads in sparse delimited file.
#' There is one version for the case in which genes are in the rows,
#' and one for the transpose case, when they are in the columns.
#' The transpose case simpler and faster.
#' For both cases we assume the first rows are the colums names and the first column are the rownames
#' @import Matrix
#' @import methods
#' @param file character(1) path to a gzipped delimited file with counts; features are rows, cells are columns
#' @param sep character(1) delimiter for fields (defaults to ',')
#' @examples
#' ## Original object
#' data("tenx_subset")
#' new <- read_sparse_csv(system.file("extdata/tenx_subset.csv.gz", package = "smallcount"))
#' identical(new, tenx_subset)
#' new <- read_sparse_csv_transpose(system.file("extdata/tenx_subset-transpose.csv.gz", package = "smallcount"))
#' identical(new, tenx_subset)
#' @export
read_sparse_csv <- function(file, sep = ","){

  conn <- gzcon(file(file.path(file), "rb"))
  ## Get cellnames
  l <- readLines(conn, n = 1)
  cns <- strsplit(l, sep)[[1]][-1]

  ## initialize lists,
  ##rns are the genenames
  ## i are the rows of non-zeros
  ## j are the columns of the non-zeros
  ## x are the non-zero values
  rns <- vector("list")
  i <- vector("list")
  j <- vector("list")
  x <- vector("list")

  ##n is counting rows
  n <- 0L
  while(TRUE){
    l <- readLines(conn, n = 1)
    if(!length(l)) break()

    n <- n + 1L

    y <- strsplit(l, ",")[[1]]
    rns[[n]] <- y[1]
    y <- as.numeric(y[-1])

    j[[n]] <- as.integer(which(y>0)) - 1L ## these are non zero columns

    i[[n]] <- rep(n, length(j[[n]])) - 1L ## these are the rows

    x[[n]] <- y[ j[[n]] + 1L] ## these are t

  }
  close(conn)

  ##turn lists into vectors
  i <- unlist(i)
  j <- unlist(j)
  x <- unlist(x)

  ##nummber of columns
  nc <- as.integer(length(cns))

  ## convert to a dgCMatrix format
  ## to do this we have to
  ## order everything in column order
  o <- order(j)
  i <- i[o]
  x <- x[o]
  ## then count for each column where it ends
  j <- c(0L, cumsum(as.integer(table(factor(j, levels = 0:(nc-1))))))

  Dim <- c(n, nc)
  Dimnames <- list(unlist(rns), cns)

  newm <- new("dgCMatrix", x = x, i = i, p = j, Dim = Dim, Dimnames = Dimnames)
  return(newm)
}


#' Read genes in columns
#' @describeIn read_sparse_csv reads csv with genes in the columns
#' @param file character(1) path to a gzipped delimited file with counts; features are rows, cells are columns
#' @param sep character(1) delimiter for fields (defaults to ',')
#' @export
read_sparse_csv_transpose <- function(file, sep = ","){

  conn <- gzcon(file(file.path(file), "rb"))
  l <- readLines(conn, n = 1)

  ## Get genenames
  rns <- strsplit(l, sep)[[1]][-1]

  ## initialize lists,
  ## cns are the cell names
  ## i is the row index of non-zeros within column
  ## x are the non-zero values
  ## p has a 0 followed by an index for each column denoting where the column ends within in i and x
  cns <- vector("list")
  i <- vector("list")
  p <- vector("list")
  x <- vector("list")

  ## n is counting rows
  n <- 0L
  ## In dgCMatrix first p is 0
  p[[1]] <- 0L

  ## Read line by line until the end
  while(TRUE){
    l <- readLines(conn, n = 1)
    if(!length(l)) break()

    n <- n + 1L

    y <- strsplit(l, sep)[[1]]
    ## this is the n-th row names
    cns[[n]] <- y[1]

    ## remove first entry which is rowname
    y <- as.numeric(y[-1])

    ## index of non zero columns
    i[[n]] <- as.integer(which(y>0)) - 1L

    ## end of the n-th column
    p[[n+1]] <- p[[n]] + length(i[[n]])

    ## non-zero values
    x[[n]] <- y[ i[[n]] + 1L] ## these are t

  }
  close(conn)

  ## Generate entries for the dgCMatrix
  i <- unlist(i)
  p <- unlist(p)
  x <- unlist(x)

  Dim <- c(as.integer(length(rns)), n)
  Dimnames <- list(rns, unlist(cns))

  ## Creat matrix, this will happen in R not in C++
  newm <- new("dgCMatrix", x = x, i = i, p = p, Dim = Dim, Dimnames = Dimnames)
  return(newm)
}

