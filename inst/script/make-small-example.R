library(Matrix)
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

fn <- "inst/extdata/tenx_subset-transpose.csv"
write.csv(t(m), file = fn, row.names = TRUE, quote = FALSE)
system(paste("gzip",fn))

tenx_subset <- as(m, "dgCMatrix")
fn <- "inst/extdata/tenx_subset.mtx"
writeMM(tenx_subset, file=fn)
system(paste("gzip",fn))


save(tenx_subset, file= "data/tenx_subset.rda", compress = "xz")

