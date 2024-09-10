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

df = pd.DataFrame(arr, index=["r1", "r2", "r3"], columns=["c1", "c2", "c3"])
csv_filename = filename + ".csv"
df.to_csv(csv_filename)
with tarfile.open(filename + ".tar.gz", "w:gz") as tar:
    tar.add(csv_filename)
os.remove(csv_filename)

target = BytesIO()
mmwrite(target, coo_matrix(arr))
mtx_filename = filename + ".mtx"
with open(mtx_filename, "wb") as f:
    f.write(target.getbuffer())
with tarfile.open(filename + ".tbz2", "w:bz2") as tar:
    tar.add(mtx_filename)
os.remove(mtx_filename)

hf = h5py.File(filename + ".h5", "w")
csc_mat = csc_matrix(arr)
hf.create_dataset("matrix/data", dtype=np.uint32, data=csc_mat.data)
hf.create_dataset("matrix/indices", dtype=np.uint32, data=csc_mat.indices)
hf.create_dataset("matrix/indptr", dtype=np.uint32, data=csc_mat.indptr)
hf.create_dataset("matrix/shape", dtype=np.uint64, data=csc_mat.shape)
hf.close()
