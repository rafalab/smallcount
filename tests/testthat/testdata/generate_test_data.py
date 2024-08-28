import h5py
from io import BytesIO
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.io import mmwrite

arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

target = BytesIO()
mmwrite(target, coo_matrix(arr))
with open("small_dense_square.mtx", "wb") as f:
    f.write(target.getbuffer())

hf = h5py.File("small_dense_square.h5", "w")
csc_mat = csc_matrix(arr)
hf.create_dataset("matrix/data", dtype=np.uint32, data=csc_mat.data)
hf.create_dataset("matrix/indices", dtype=np.uint32, data=csc_mat.indices)
hf.create_dataset("matrix/indptr", dtype=np.uint32, data=csc_mat.indptr)
hf.create_dataset("matrix/shape", dtype=np.uint64, data=csc_mat.shape)
hf.close()
