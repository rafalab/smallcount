#ifndef SMALLCOUNT_MTX_FILE_READER_H_
#define SMALLCOUNT_MTX_FILE_READER_H_

#include <fstream>

#include "sparse_matrix.h"
#include "tenx_file_params.h"

namespace smallcount {

// File reader to construct sparse matrices from .mtx files.
class MtxFileReader {
   public:
    // Reads the contents of an .mtx file into a SparseMatrix object, labelling
    // the rows and columns with features and barcodes, respectively.
    static void read(std::ifstream &matrix_file, std::ifstream &barcodes_file,
                     std::ifstream &features_file, const TenxFileParams &params,
                     SparseMatrix &matrix);

   private:
    // Static class. Should not be instantiated.
    MtxFileReader() = default;
};

}  // namespace smallcount

#endif
