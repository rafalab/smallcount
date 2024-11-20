"""Generates test data for the read_sparse_matrix function.

Test matrix:
   r1 r2 r3
c1  1  2  3
c2  4  5  6
c3  7  8  9

Output file formats:
    1. .csv, tarballed and gzipped
    2. .mtx, tarballed and bzipped
    3. .h5
"""

import h5py
from io import BytesIO
import numpy as np
import os
import pandas as pd
from scipy.sparse import coo_matrix, csc_matrix
from scipy.io import mmwrite
import tarfile

filename = "small_dense_square"
arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

# .csv file
rows = ["r1", "r2", "r3"]
cols = ["c1", "c2", "c3"]
df = pd.DataFrame(arr, index=rows, columns=cols)
csv_filename = filename + ".csv"
df.to_csv(csv_filename)
with tarfile.open(filename + ".tar.gz", "w:gz") as tar:
    tar.add(csv_filename)
os.remove(csv_filename)

# .mtx file
target = BytesIO()
mmwrite(target, coo_matrix(arr))
mtx_filename = filename + ".mtx"
with open(mtx_filename, "wb") as f:
    f.write(target.getbuffer())
with tarfile.open(filename + ".tbz2", "w:bz2") as tar:
    tar.add(mtx_filename)
os.remove(mtx_filename)

# .h5 file
hf = h5py.File(filename + ".h5", "w")
csc_mat = csc_matrix(arr)
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

hf.create_dataset("genome/data", dtype=np.uint32, data=csc_mat.data)
hf.create_dataset("genome/indices", dtype=np.uint32, data=csc_mat.indices)
hf.create_dataset("genome/indptr", dtype=np.uint32, data=csc_mat.indptr)
hf.create_dataset("genome/shape", dtype=np.uint64, data=csc_mat.shape)
hf.create_dataset("genome/gene_names", data=utf8_rows)
hf.create_dataset("genome/barcodes", data=utf8_cols)
hf.close()
