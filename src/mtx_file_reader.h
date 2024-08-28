#ifndef SMALLCOUNT_MTX_FILE_READER_H_
#define SMALLCOUNT_MTX_FILE_READER_H_

#include <fstream>

#include "sparse_matrix.h"

namespace smallcount {

// File reader to construct sparse matrices from .mtx files.
class MtxFileReader {
   public:
    // Reads the contents of an .mtx file into a SparseMatrix object.
    static void read(std::ifstream &file, SparseMatrix &matrix);

   private:
    // Static class. Should not be instantiated externally.
    MtxFileReader() = default;
};

}  // namespace smallcount

#endif
