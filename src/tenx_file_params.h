#ifndef SMALLCOUNT_TENX_FILE_PARAMS_H_
#define SMALLCOUNT_TENX_FILE_PARAMS_H_

#include <optional>
#include <string>

namespace smallcount {

// Parameters specifying how to read the CellRanger output files for 10x
// Genomics data.
struct TenxFileParams {
    // Whether to use the cell barcodes as column names.
    bool use_barcode_col_names;
    // Whether to use the Ensembl IDs as row names. Otherwise use gene symbols.
    bool use_id_row_names;

    // FOR MTX FILES:
    // Whether to use features.tsv for features. Otherwise use genes.tsv.
    bool use_features_tsv;

    // FOR HDF5 FILES:
    // Name of the HDF5 group containing the matrix datasets for CellRanger v2.
    std::optional<std::string> genome = std::nullopt;
};

}  // namespace smallcount

#endif
