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
pc <- pca_poisson_residuals(tenx_subset)
```

This will give same results as 
```
x <- t(as.matrix(tenx_subset))
pc_old <- prcomp(x)
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
