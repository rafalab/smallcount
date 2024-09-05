import h5py
from io import BytesIO
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csc_matrix
from scipy.io import mmwrite

filename = "small_dense_square"
arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

df = pd.DataFrame(arr, index=["r1", "r2", "r3"], columns=["c1", "c2", "c3"])
df.to_csv(filename + ".csv")

target = BytesIO()
mmwrite(target, coo_matrix(arr))
with open(filename + ".mtx", "wb") as f:
    f.write(target.getbuffer())

hf = h5py.File(filename + ".h5", "w")
csc_mat = csc_matrix(arr)
hf.create_dataset("matrix/data", dtype=np.uint32, data=csc_mat.data)
hf.create_dataset("matrix/indices", dtype=np.uint32, data=csc_mat.indices)
hf.create_dataset("matrix/indptr", dtype=np.uint32, data=csc_mat.indptr)
hf.create_dataset("matrix/shape", dtype=np.uint64, data=csc_mat.shape)
hf.close()
