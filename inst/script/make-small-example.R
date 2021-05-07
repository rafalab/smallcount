library(Matrix)
library(TENxBrainData)
tenx <- TENxBrainData()

set.seed(2020-05-07)
rind <- sample(nrow(tenx), 100)
cind <- sample(nrow(tenx), 10000)

m <- as.matrix(counts(tenx[rind,cind]))
rownames(m) <- rowData(tenx)$Ensembl[rind]
fn <- "inst/extdata/small_example.csv"
write.csv(m, file = fn, row.names = TRUE, quote = FALSE)
system(paste("gzip",fn))

fn <- "inst/extdata/small_example-transpose.csv"
write.csv(t(m), file = fn, row.names = TRUE, quote = FALSE)
system(paste("gzip",fn))

small_example <- as(m, "dgCMatrix")
fn <- "inst/extdata/small_example.mtx"
writeMM(small_example, file="inst/extdata/small_example.mtx")
system(paste("gzip",fn))


saveRDS(small_example, file= "inst/extdata/small_example.rds")

