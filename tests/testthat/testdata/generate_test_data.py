"""Generates test data for the read_sparse_matrix function.

Test matrix:
   r1 r2 r3
c1  1  2  3
c2  4  5  6
c3  7  8  9

Output file formats:
    1. .csv, tarballed and gzipped
    2. .mtx, tarballed and bzipped
    3. .h5, using CellRanger v3 format
    4. .h5, using CellRanger v2 format
"""

import h5py
from io import BytesIO
import numpy as np
import os
import pandas as pd
from scipy.sparse import coo_matrix, csc_matrix
from scipy.io import mmwrite
import tarfile

filedir = "tests/testthat/testdata/"
filename = "small_dense_square"
arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

# .csv file
rows = ["r1", "r2", "r3"]
cols = ["c1", "c2", "c3"]
df = pd.DataFrame(arr, index=rows, columns=cols)
csv_filename = filename + ".csv"
df.to_csv(csv_filename)
with tarfile.open(filedir + filename + ".tar.gz", "w:gz") as tar:
    tar.add(csv_filename)
os.remove(csv_filename)


# .mtx directory (Cell Ranger v3)
target = BytesIO()
mmwrite(target, coo_matrix(arr))
tempdir = filename + "_v3/"
os.mkdir(tempdir)
mtx_filename = tempdir + "matrix.mtx"
with open(mtx_filename, "wb") as f:
    f.write(target.getbuffer())
barcodes_filename = tempdir + "barcodes.tsv"
pd.DataFrame(cols).to_csv(barcodes_filename, sep="\t", header=False,
                          index=False)
features_filename = tempdir + "features.tsv"
# Put row information in IDs.
features = pd.DataFrame(np.array([rows, cols]).T)
features.to_csv(features_filename, sep="\t", header=False, index=False)
with tarfile.open(filedir + tempdir[:-1] + ".tbz2", "w:bz2") as tar:
    tar.add(mtx_filename)
    tar.add(barcodes_filename)
    tar.add(features_filename)
os.remove(mtx_filename)
os.remove(barcodes_filename)
os.remove(features_filename)
os.rmdir(tempdir)

# .mtx directory (Cell Ranger v2)
prefix = "prefix_"
prefix_dir = filedir + filename + "_v2/"
if not os.path.exists(prefix_dir):
    os.mkdir(prefix_dir)
mtx_filename2 = prefix_dir + prefix + "matrix.mtx"
with open(mtx_filename2, "wb") as f:
    f.write(target.getbuffer())
barcodes_filename2 = prefix_dir + prefix + "barcodes.tsv"
pd.DataFrame(cols).to_csv(barcodes_filename2, sep="\t", header=False,
                          index=False)
features_filename2 = prefix_dir + prefix + "genes.tsv"
# Put row information in names/symbols.
features2 = pd.DataFrame(np.array([cols, rows, cols]).T)
features2.to_csv(features_filename2, sep="\t", header=False, index=False)


# .h5 file (Cell Ranger v3)
csc_mat = csc_matrix(arr)

hf = h5py.File(filedir + filename + "_v3.h5", "w")
hf.create_dataset("matrix/data", dtype=np.uint32, data=csc_mat.data)
hf.create_dataset("matrix/indices", dtype=np.uint32, data=csc_mat.indices)
hf.create_dataset("matrix/indptr", dtype=np.uint32, data=csc_mat.indptr)
hf.create_dataset("matrix/shape", dtype=np.uint64, data=csc_mat.shape)

str_len = 10
to_utf8 = lambda x: x.encode("utf-8")
utf8_rows = np.array([to_utf8(r) for r in rows], dtype=f"S{str_len}")
utf8_cols = np.array([to_utf8(c) for c in cols], dtype=f"S{str_len}")
hf.create_dataset("matrix/features/id", data=utf8_rows)
hf.create_dataset("matrix/barcodes", data=utf8_cols)
hf.close()

# .h5 file (Cell Ranger v2)
hf2 = h5py.File(filedir + filename + "_v2.h5", "w")
hf2.create_dataset("genome/data", dtype=np.uint32, data=csc_mat.data)
hf2.create_dataset("genome/indices", dtype=np.uint32, data=csc_mat.indices)
hf2.create_dataset("genome/indptr", dtype=np.uint32, data=csc_mat.indptr)
hf2.create_dataset("genome/shape", dtype=np.uint64, data=csc_mat.shape)
hf2.create_dataset("genome/gene_names", data=utf8_rows)
hf2.create_dataset("genome/barcodes", data=utf8_cols)
hf2.close()
