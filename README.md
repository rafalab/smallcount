## R package with Methods for Small Counts Stored in a Sparse Matrix 

You can install like this:

```
devtools::install_github("rafalab/smallcount")
```

Make sure your matrix is of class `dgCMatrix`

```
library(smallcount)
data(tenx_subset)
class(tenx_subset)
```

You can compute Pearson residuals followed by PCA like this:

```
dim(tenx_subset)
system.time({pc <- pca_poisson_residuals(tenx_subset, residual = "pearson")})
```

This will give same results as 
```
safe_divide <- function(a, b) {
    ifelse(b == 0, 0, a / b)
}
system.time({
    n <- colSums(tenx_subset)
    rate <- rowSums(tenx_subset) / sum(n)
    rate_n <- outer(rate, n)
    pearson_residuals <- safe_divide(tenx_subset - rate_n, sqrt(rate_n))
    pc_old <- prcomp(t(pearson_residuals), center = FALSE, rank. = 50)
})
```

but about 50 times faster and using less memory.


## Requirements

The matrix should store raw counts not transformed values. Genes should be in the rows and cells in the columns. If stored as regular matrix, convert to sparse one like this:

```
x <- matrix(rpois(1000000,0.1), 100, 1000)
y <- as(x, "dgCMatrix")
```

You can see the gains in memory consumption like this:

```
object.size(x)
object.size(y)
```
