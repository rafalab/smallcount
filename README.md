## R package with Methods for Small Counts Stored in a Sparse Matrix 

You can install like this:

```
devtools::install_github("rafalab/smallcount")
```

Make sure your matrix is of class `SVT_SparseMatrix`:

```
library(smallcount)
data(tenx_subset)
class(tenx_subset)
```

You can compute Pearson residuals followed by PCA like this:

```
dim(tenx_subset)
system.time({pc <- poissonPca(tenx_subset, transform = "pearson")})
```

This will give same results as:
```
x <- as.matrix(tenx_subset)
safe_divide <- function(a, b) {
    ifelse(b == 0, 0, a / b)
}
system.time({
    mu_hat <- outer(rowSums(x) / sum(x), colSums(x))
    pearson_residuals <- safe_divide(x - mu_hat, sqrt(mu_hat))
    pc_old <- prcomp(t(pearson_residuals), center = FALSE, rank. = 50)
})
```

but about 50 times faster and using less memory.


## Requirements

The matrix should store raw counts, not transformed values. Genes should be in the rows and cells in the columns. If stored as a regular matrix, convert to a sparse one like this:

```
x <- matrix(rpois(100000, 0.1), 100, 1000)
y <- as(x, "SparseMatrix")
```

You can see the gains in memory consumption like this:

```
object.size(x)
object.size(y)
```
