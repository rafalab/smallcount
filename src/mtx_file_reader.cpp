#include "mtx_file_reader.h"

#include <cstdlib>
#include <fstream>
#include <optional>
#include <string>

#include "Rcpp.h"
#include "sparse_matrix.h"
#include "tenx_file_params.h"

using namespace Rcpp;

namespace smallcount {
namespace {

// Non-zero entry in a sparse matrix.
struct MatrixData {
    int row;  // Row index (1-based index, for consistency with R)
    int col;  // Column index (1-based index, for consistency with R)
    int val;  // Value
};

// Line of data/metadata in an .mtx file.
// Typically corresponds to a non-zero entry in a sparse matrix.
// Contains matrix metadata for the first line of the .mtx file.
struct MtxLine {
    // Entry row index (or total number of rows).
    int row;
    // Entry column index (or total number of columns).
    int col;
    // Entry value (or total number of non-zero values).
    size_t val;

    MatrixData data() const {
        return {.row = row, .col = col, .val = static_cast<int>(val)};
    }
    MatrixMetadata metadata() const {
        return {.nrow = row, .ncol = col, .nval = val};
    }
};

// Reads the first or second column from a .tsv file.
std::vector<std::string> readTsvNames(std::ifstream &file,
                                      const std::string &name, int size,
                                      bool first_row) {
    std::string line;
    std::vector<std::string> names;
    names.reserve(size);
    while (std::getline(file, line)) {
        size_t tab_index = line.find('\t');
        if (tab_index == std::string::npos && !first_row) {
            stop("Invalid features/genes. Could not locate second column.");
        }
        size_t start_index = first_row ? 0 : tab_index + 1;
        size_t end_index =
            first_row ? tab_index : line.find('\t', tab_index + 1);
        names.emplace_back(line.substr(start_index, end_index - start_index));
    }
    if (names.size() != size) {
        warning(
            "Number of %s does not match the specifications in the metadata "
            "(%zu vs. %d)",
            name, names.size(), size);
    }
    return names;
}

// Generates metadata with row and column names.
MatrixMetadata createMetadata(MtxLine entry, std::ifstream &barcodes_file,
                              std::ifstream &features_file,
                              const TenxFileParams &params) {
    MatrixMetadata metadata = entry.metadata();
    auto row_names =
        readTsvNames(features_file, "features/genes", metadata.nrow,
                     /*first_row=*/params.use_id_row_names);
    std::vector<std::string> col_names;
    if (params.use_barcode_col_names) {
        col_names = readTsvNames(barcodes_file, "barcodes", metadata.ncol,
                                 /*first_row=*/true);
    }
    metadata.row_names = std::move(row_names);
    metadata.col_names = std::move(col_names);
    return metadata;
}

// Parses a single line of an .mtx file.
std::optional<MtxLine> parseMtxLine(const std::string &line, size_t line_num) {
    // Ignore comments.
    if (line[0] == '%') {
        return std::nullopt;
    }
    // Read row, column, and value information.
    MtxLine result;
    const char *c_str = line.c_str();
    char *str_end;
    result.row = strtol(c_str, &str_end, /*__base=*/10);
    result.col = strtol(str_end, &str_end, /*__base=*/10);
    result.val = strtol(str_end, &str_end, /*__base=*/10);
    if (result.row == 0 || result.col == 0 || result.val == 0) {
        stop(
            "Unexpected entry. Line %zu does not specify three positive "
            "integers:\n%s",
            line_num, line);
        return std::nullopt;
    }
    return result;
}

}  // namespace

SvtSparseMatrix MtxFileReader::read(std::ifstream &matrix_file,
                                    std::ifstream &barcodes_file,
                                    std::ifstream &features_file,
                                    const TenxFileParams &params) {
    Svt svt;
    MatrixMetadata metadata;

    std::string line;
    size_t line_num = 0;
    size_t non_zero_count = 0;
    while (std::getline(matrix_file, line)) {
        line_num++;
        const auto entry = parseMtxLine(line, line_num);
        if (!entry.has_value()) {
            continue;
        } else if (svt.empty()) {
            metadata =
                createMetadata(*entry, barcodes_file, features_file, params);
            svt = Svt(metadata.ncol, SvtEntry(2));
        } else {
            svt[entry->col - 1][kSvtRowInd].emplace_back(entry->row - 1);
            svt[entry->col - 1][kSvtValInd].emplace_back(entry->val);
            non_zero_count++;
        }
    }

    if (non_zero_count != metadata.nval) {
        stop(
            "Inconsistent entry count. Number of non-zero entries does not "
            "match the total specified in the matrix metadata (%zu != %zu).",
            non_zero_count, metadata.nval);
    }

    return SvtSparseMatrix(std::move(svt), std::move(metadata));
}

}  // namespace smallcount
