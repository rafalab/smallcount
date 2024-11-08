#include "mtx_file_reader.h"

#include <cstdlib>
#include <fstream>
#include <optional>
#include <string>

#include "Rcpp.h"
#include "sparse_matrix.h"

using namespace Rcpp;

namespace smallcount {
namespace {

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

void MtxFileReader::read(std::ifstream &file, SparseMatrix &matrix) {
    std::string line;
    size_t line_num = 0;
    size_t non_zero_count = 0;
    bool is_initialized = false;
    while (std::getline(file, line)) {
        line_num++;
        const auto entry = parseMtxLine(line, line_num);
        if (!entry.has_value()) {
            continue;
        } else if (is_initialized) {
            matrix.addEntry(entry->data());
            non_zero_count++;
        } else {
            matrix.init(entry->metadata());
            is_initialized = true;
        }
    }
    if (non_zero_count != matrix.nval()) {
        stop(
            "Inconsistent entry count. Number of non-zero entries does not "
            "match the total specified in the matrix metadata (%zu != %zu).",
            non_zero_count, matrix.nval());
    }
}

}  // namespace smallcount
