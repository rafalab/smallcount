#include "mtx_file_reader.h"

#include <cstdlib>
#include <fstream>
#include <optional>
#include <string>

#include "Rcpp.h"
#include "mtx_file.h"
#include "sparse_matrix.h"

using namespace Rcpp;

namespace smallcount {
namespace {

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
        warning(
            "Error processing line %zu. Encountered zeroes or misconfigured "
            "entry: %d %d %d",
            line_num, result.row, result.col, result.val);
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
        }
        if (is_initialized) {
            matrix.addEntry(entry->data());
            non_zero_count++;
        } else {
            matrix.init(entry->metadata());
            is_initialized = true;
        }
    }
    if (non_zero_count != matrix.nval()) {
        warning(
            "Count of non-zero data does not match the count specified in "
            "the matrix metadata (%zu != %zu).",
            non_zero_count, matrix.nval());
    }
}

}  // namespace smallcount
