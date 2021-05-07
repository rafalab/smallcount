### If the genes are represented by rows
read_sparse_csv <- function(fn){
  require(Matrix)
  conn <- gzcon(file(file.path(fn), "rb"))
  l <- readLines(conn, n = 1)
  cns <- strsplit(l, ",")[[1]][-1]
  rns <- vector("list")
  i <- vector("list")
  j <- vector("list")
  x <- vector("list")
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

  i <- unlist(i)
  j <- unlist(j)
  x <- unlist(x)

  nc <- as.integer(length(cns))
  o <- order(j)
  i <- i[o]
  x <- x[o]
  j <- c(0L, cumsum(as.integer(table(factor(j, levels = 0:(nc-1))))))

  Dim <- c(n, nc)
  Dimnames <- list(unlist(rns), cns)

  newm <- new("dgCMatrix", x = x, i = i, p = j, Dim = Dim, Dimnames = Dimnames)
  return(newm)
}



### If the genes are represented by columns
read_sparse_csv_transpose <- function(fn){
  require(Matrix)
  conn <- gzcon(file(file.path(fn), "rb"))
  l <- readLines(conn, n = 1)
  rns <- strsplit(l, ",")[[1]][-1]
  cns <- vector("list")
  i <- vector("list")
  p <- vector("list")
  x <- vector("list")
  n <- 0L
  p[[1]] <- 0L
  while(TRUE){
    l <- readLines(conn, n = 1)
    if(!length(l)) break()

    n <- n + 1L

    y <- strsplit(l, ",")[[1]]
    cns[[n]] <- y[1]
    y <- as.numeric(y[-1])

    i[[n]] <- as.integer(which(y>0)) - 1L ## these are non zero rows

    p[[n+1]] <- p[[n]] + length(i[[n]])

    x[[n]] <- y[ i[[n]] + 1L] ## these are t

  }
  close(conn)

  i <- unlist(i)
  p <- unlist(p)
  x <- unlist(x)

  Dim <- c(as.integer(length(rns)), n)
  Dimnames <- list(rns, unlist(cns))

  newm <- new("dgCMatrix", x = x, i = i, p = p, Dim = Dim, Dimnames = Dimnames)
  return(newm)
}


## Test
m <- readRDS("inst/extdata/small_example.rds")

newm <- read_sparse_csv("inst/extdata/small_example.csv.gz")
identical(newm, m)

newm <- read_sparse_csv_transpose("inst/extdata/small_example-transpose.csv.gz")
identical(newm, m)


