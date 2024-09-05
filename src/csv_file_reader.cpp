#include "csv_file_reader.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "Rcpp.h"
#include "sparse_matrix.h"

using namespace Rcpp;

namespace smallcount {
namespace {

struct NonZeroEntry {
    int col;  // Column index (1-based)
    int val;  // Value
};

using NonZeroEntries = std::vector<NonZeroEntry>;

// Returns the non-zero entries in a line of .csv file.
NonZeroEntries parseCsvLine(const std::string &line, int ncol, int line_num) {
    // Ignore the row name.
    const char *val_start = strchr(line.c_str(), ',');
    if (val_start == nullptr) {
        stop("No comma found on line %d.", line_num);
    }
    char *val_end;
    NonZeroEntries entries;
    int col = 0;
    while (*val_start != '\0') {
        col++;
        // Read the next value, skipping over the comma.
        int val = strtol(val_start + 1, &val_end, /*__base=*/10);
        val_start = val_end;
        if (val != 0) {
            entries.emplace_back(NonZeroEntry{.col = col, .val = val});
        }
    }
    if (col != ncol) {
        warning(
            "Error processing line %d. Expected %d entries but encountered %d.",
            line_num, ncol, col);
    }
    return entries;
}

}  // namespace

void CsvFileReader::read(std::ifstream &file, SparseMatrix &matrix) {
    std::string line;
    int line_num = 0;
    int ncol = 0;
    size_t nval = 0;
    std::vector<NonZeroEntries> sparse_mat;
    while (std::getline(file, line)) {
        line_num++;
        if (line_num > 1) {
            auto nz_entries = parseCsvLine(line, ncol, line_num);
            nval += nz_entries.size();
            sparse_mat.emplace_back(std::move(nz_entries));
        } else {
            // When row names are provided, the number of columns is given by
            // the number of commas.
            ncol = std::count(line.begin(), line.end(), ',');
        }
    }
    matrix.init(
        MatrixMetadata{.nrow = line_num - 1, .ncol = ncol, .nval = nval});
    for (int row = 0; row < sparse_mat.size(); row++) {
        for (const NonZeroEntry &entry : sparse_mat[row]) {
            // Shift row by +1 for 1-based indexing.
            matrix.addEntry(
                MatrixData{.row = row + 1, .col = entry.col, .val = entry.val});
        }
    }
}

}  // namespace smallcount
