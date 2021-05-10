group_means <- function(y, g){

  if(!is.factor(g)){
    warning("Coercing g into a factor")
    g <- as.factor(g)
  }

  js <- as.numeric(g)

  n <- sapply(split(seq_along(g), g), length)

  m <- matrix(0, nrow(y),  length(n))

  for(j in 1:ncol(y)){
    ind <- (y@p[j]+1):y@p[j+1]
    real_ind <- y@i[ind] + 1
    k <- js[j]
    m[real_ind, k] <- m[real_ind, k] + y@x[ind]/n[k]
  }
  colnames(m) <- levels(g)
  rownames(m) <- rownames(y)

  return(m)
}

