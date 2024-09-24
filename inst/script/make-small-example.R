library(Matrix)
library(SparseArray)
library(TENxBrainData)
tenx <- TENxBrainData()

set.seed(2020-05-08)
rind <- sort(sample(nrow(tenx), 500))
cind <- sort(sample(nrow(tenx), 10000))

m <- as.matrix(counts(tenx[rind,cind]))
rownames(m) <- rowData(tenx)$Ensembl[rind]
fn <- "inst/extdata/tenx_subset.csv"
write.csv(m, file = fn, row.names = TRUE, quote = FALSE)
system(paste("gzip",fn))

tenx_subset_dgc <- as(m, "dgCMatrix")
fn <- "inst/extdata/tenx_subset.mtx"
writeMM(tenx_subset_dgc, file=fn)
system(paste("gzip",fn))

tenx_subset <- as(m, "SparseMatrix")
save(tenx_subset, file= "data/tenx_subset.rda", compress = "xz")
