#ifndef SMALLCOUNT_CSV_FILE_READER_H_
#define SMALLCOUNT_CSV_FILE_READER_H_

#include <fstream>

#include "sparse_matrix.h"

namespace smallcount {

// File reader to construct sparse matrices from .csv files.
class CsvFileReader {
   public:
    // Converts the contents of a .csv file into an SvtSparseMatrix.
    static SvtSparseMatrix read(std::ifstream &file);

   private:
    // Static class. Should not be instantiated.
    CsvFileReader() = default;
};

}  // namespace smallcount

#endif
