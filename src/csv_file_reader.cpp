#include "csv_file_reader.h"

#include <algorithm>
#include <cmath>
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
        stop("No comma delimiter on line 1.");
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

// Reads the row name and non-zero data given in a line of a .csv file,
// directly adding entries to the provided SVT matrix.
void parseCsvLine(const std::string &line, int ncol, int line_num, int row_idx,
                  Svt &svt, std::string &row_name, size_t *nval) {
    // Read the row name.
    const char *val_start = strchr(line.c_str(), ',');
    if (val_start == nullptr) {
        stop("No comma delimiter on line %d.", line_num);
    }
    row_name = std::string(line.c_str(), val_start);

    // Read the row data.
    char *val_end;
    int col = 0;
    while (*val_start != '\0') {
        // Read the next value, skipping over the comma.
        double val = strtof(val_start + 1, &val_end);
        while (*val_end != '\0' && *val_end != ',') {
            val_end++;
        }
        if (trunc(val) != val) {
            stop(
                "Unexpected float value. Expected integer in row %d, column %d "
                "but encountered a floating-point number: %f.",
                line_num, col + 2, val);
        } else if (val != 0) {
            svt[col][kSvtRowInd].emplace_back(row_idx);
            svt[col][kSvtValInd].emplace_back(val);
            (*nval)++;
        }
        val_start = val_end;
        col++;
    }
    if (col != ncol) {
        stop(
            "Inconsistent column count. Expected %d columns (from header) but "
            "encountered %d in row %d.",
            ncol, col, line_num);
    }
}

}  // namespace

SvtSparseMatrix CsvFileReader::read(std::ifstream &file) {
    Svt svt;
    std::vector<std::string> col_names;
    std::vector<std::string> row_names;
    int ncol = 0;
    size_t nval = 0;

    std::string line;
    int line_num = 0;
    while (std::getline(file, line)) {
        line_num++;
        if (line_num > 1) {
            // Initialize each column with two vectors (rows and values)
            if (svt.empty()) svt = Svt(ncol, SvtEntry(2));
            std::string row_name;
            parseCsvLine(line, /*ncol=*/ncol, /*line_num=*/line_num,
                         /*row_idx=*/line_num - 2, svt, row_name, &nval);
            row_names.emplace_back(row_name);
        } else {
            col_names = readColumnNames(line);
            ncol = col_names.size();
        }
    }

    MatrixMetadata metadata{.nrow = line_num - 1,
                            .ncol = ncol,
                            .nval = nval,
                            .row_names = std::move(row_names),
                            .col_names = std::move(col_names)};
    return SvtSparseMatrix(std::move(svt), std::move(metadata));
}

}  // namespace smallcount
