## R package with Methods for Small Counts Stored in a Sparse Matrix 

You can instlal like this:

```
devtools::install_github("rafalab/smallcount")
```

Make sure your matrix is of class `dgCMatrix`

```
library(smallcount)
data(tenx_subset)
class(tenx_subset)
```

Compute Pearson residuals followed by PCA:

```
dim(tenx_subset)
system.time({pc <- pca_poisson_residuals(tenx_subset)})
```

This will give same results as 
```
x <- t(as.matrix(tenx_subset))
system.time({
pc_old <- prcomp(x)
})

```

but about 50 times faster and using less memory.


## Requirements

The matrix should store raw counts not transoformed values. Genes should be in the rows and cells in the columns. If stored as regular matrix, convert to sparse one like this:

```
x <- matrix(rpois(1000000,0.1), 100, 1000)
object.size(x)
x <- as(x, "dgCMatrix")
object.size(x)
```

