#ifndef SMALLCOUNT_FILE_PARAMS_H_
#define SMALLCOUNT_FILE_PARAMS_H_

#include <optional>
#include <string>

namespace smallcount {

struct FileParams {
    // Whether to use the cell barcodes as column names.
    bool use_barcode_col_names;
    // Whether to use the Ensembl IDs as row names. Otherwise use gene symbols.
    bool use_id_row_names;
    // Name of the HDF5 group containing the matrix datasets for CellRanger v2.
    std::optional<std::string> genome = std::nullopt;
};

}  // namespace smallcount

#endif
