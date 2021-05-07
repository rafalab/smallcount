### The functions are prototypes for a parser of csv files
### that returns the components needed to form  dgCMatrix
### There is one version for the case in which genes are in the columns of the csv
### and one for when they are in the rows. The former is simpler and faster so I include it first.
### for both we assume the first rows are the colums names and the first column are the rownames

### Genes are represented by columns
read_sparse_csv_transpose <- function(fn){
  require(Matrix)

  conn <- gzcon(file(file.path(fn), "rb"))
  l <- readLines(conn, n = 1)

  ## Get genenames
  rns <- strsplit(l, ",")[[1]][-1]


  ## initialize lists,
  ##cns are the cell names
  ## i is the columns index of non-zeros within row,
  ## p is how many nonzero in each column
  ## x are the non-zero values
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

    y <- strsplit(l, ",")[[1]]
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

### Genes in rows
read_sparse_csv <- function(fn){
  require(Matrix)

  conn <- gzcon(file(file.path(fn), "rb"))
  ## Get cellnames
  l <- readLines(conn, n = 1)
  cns <- strsplit(l, ",")[[1]][-1]

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
  ## then count for each column how many entries in each columns
  o <- order(j)
  i <- i[o]
  x <- x[o]
  j <- c(0L, cumsum(as.integer(table(factor(j, levels = 0:(nc-1))))))

  Dim <- c(n, nc)
  Dimnames <- list(unlist(rns), cns)

  newm <- new("dgCMatrix", x = x, i = i, p = j, Dim = Dim, Dimnames = Dimnames)
  return(newm)
}




## Test
m <- readRDS("inst/extdata/small_example.rds")

newm <- read_sparse_csv("inst/extdata/small_example.csv.gz")
identical(newm, m)

newm <- read_sparse_csv_transpose("inst/extdata/small_example-transpose.csv.gz")
identical(newm, m)


