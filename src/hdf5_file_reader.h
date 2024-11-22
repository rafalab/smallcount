#ifndef SMALLCOUNT_HDF5_FILE_READER_H_
#define SMALLCOUNT_HDF5_FILE_READER_H_

#include <optional>
#include <string>

#include "hdf5.h"
#include "sparse_matrix.h"
#include "tenx_file_params.h"

namespace smallcount {

// File reader to construct sparse matrices from Cell Ranger HDF5 files.
class Hdf5FileReader {
   public:
    // Reads the contents of an HDF5 file into a SparseMatrix object.
    static void read(hid_t file, const TenxFileParams &params,
                     SparseMatrix &matrix);

   private:
    // Static class. Should not be instantiated.
    Hdf5FileReader() = default;
};

}  // namespace smallcount

#endif
