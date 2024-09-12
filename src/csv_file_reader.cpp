#include "csv_file_reader.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "Rcpp.h"
#include "sparse_matrix.h"

using namespace Rcpp;

namespace smallcount {
namespace {

// Non-zero entry in a row of a .csv file.
struct NonZeroEntry {
    int col;  // Column index (1-based)
    int val;  // Value
};

using NonZeroEntries = std::vector<NonZeroEntry>;

// Row of data in a .csv file.
struct CsvLine {
    std::string row_name;
    NonZeroEntries nz_entries;
};

// Returns the column names given in the first line of a .csv file.
std::vector<std::string> readColumnNames(const std::string &line) {
    std::vector<std::string> col_names;
    // Skip the top-left corner of the .csv file.
    const char *val_start = strchr(line.c_str(), ',');
    if (val_start == nullptr) {
        stop("No comma found on line 1.");
    }
    // Read the column names.
    while (true) {
        val_start++;
        const char *val_end = strchr(val_start, ',');
        if (val_end != nullptr) {
            col_names.emplace_back(val_start, val_end);
        } else {
            col_names.emplace_back(val_start);
            break;
        }
        val_start = val_end;
    }
    return col_names;
}

// Reads the row name and non-zero data given in a line of a .csv file.
CsvLine parseCsvLine(const std::string &line, int ncol, int line_num) {
    // Read the row name.
    CsvLine row;
    const char *val_start = strchr(line.c_str(), ',');
    if (val_start == nullptr) {
        stop("No comma found on line %d.", line_num);
    }
    row.row_name = std::string(line.c_str(), val_start);

    // Read the row data.
    char *val_end;
    int col = 0;
    while (*val_start != '\0') {
        col++;
        // Read the next value, skipping over the comma.
        int val = strtol(val_start + 1, &val_end, /*__base=*/10);
        val_start = val_end;
        if (val != 0) {
            row.nz_entries.emplace_back(NonZeroEntry{.col = col, .val = val});
        }
    }
    if (col != ncol) {
        warning(
            "Error processing line %d. Expected %d entries but encountered %d.",
            line_num, ncol, col);
    }
    return row;
}

}  // namespace

void CsvFileReader::read(std::ifstream &file, SparseMatrix &matrix) {
    std::string line;
    int line_num = 0;
    MatrixMetadata metadata{.nval = 0};
    std::vector<NonZeroEntries> sparse_mat;
    while (std::getline(file, line)) {
        line_num++;
        if (line_num > 1) {
            auto row = parseCsvLine(line, metadata.ncol, line_num);
            metadata.nval += row.nz_entries.size();
            metadata.row_names.push_back(row.row_name);
            sparse_mat.emplace_back(std::move(row.nz_entries));
        } else {
            metadata.col_names = readColumnNames(line);
            metadata.ncol = metadata.col_names.size();
        }
    }
    metadata.nrow = line_num - 1;

    matrix.init(std::move(metadata));
    for (int row = 0; row < sparse_mat.size(); row++) {
        for (const NonZeroEntry &entry : sparse_mat[row]) {
            // Shift row by +1 for 1-based indexing.
            matrix.addEntry(
                MatrixData{.row = row + 1, .col = entry.col, .val = entry.val});
        }
    }
}

}  // namespace smallcount
