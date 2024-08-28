#include "hdf5_file_reader.h"

#include <cstdint>
#include <string>
#include <vector>

#include "Rcpp.h"
#include "hdf5.h"
#include "sparse_matrix.h"

using namespace Rcpp;

namespace smallcount {
namespace {

static constexpr char kNonZeroRows[] = "matrix/indices";
static constexpr char kNonZeroData[] = "matrix/data";
static constexpr char kColumnIndices[] = "matrix/indptr";

static constexpr char kDims[] = "matrix/shape";

// Helper function to map C++ types to HDF5 PredTypes.
template <typename T>
hid_t h5Type();
template <>
hid_t h5Type<uint32_t>() {
    return H5T_NATIVE_UINT32;
}
template <>
hid_t h5Type<uint64_t>() {
    return H5T_NATIVE_UINT64;
}

// Reads the contents of an HDF5 dataset into a vector of type T.
template <typename T>
std::vector<T> readDataset(hid_t file, const std::string &dataset_name) {
    hid_t dataset = H5Dopen2(file, dataset_name.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t num_entries;
    H5Sget_simple_extent_dims(dataspace, &num_entries, NULL);
    std::vector<T> data(num_entries);
    H5Dread(dataset, h5Type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

    H5Sclose(dataspace);
    H5Dclose(dataset);
    return data;
}

// Returns the number of rows and columns in the matrix, respectively.
std::vector<uint64_t> getDimensions(hid_t file) {
    std::vector<uint64_t> dims = readDataset<uint64_t>(file, kDims);
    if (dims.size() != 2) {
        stop("Dataset \"%s\" has %zu entries instead of 2.", kDims,
             dims.size());
    }
    return dims;
}

// Checks if a matrix entry is valid for a given matrix.
void checkValid(const SparseMatrix &matrix, const MatrixData &entry) {
    const int row = entry.row;
    const int col = entry.col;
    if (row <= 0 || row > matrix.nrow() || col <= 0 || col > matrix.ncol()) {
        stop("Coordinate (%d, %d) is out of bounds for a %d x %d matrix", row,
             col, matrix.nrow(), matrix.ncol());
    }
}

}  // namespace

void Hdf5FileReader::read(hid_t file, SparseMatrix &matrix) {
    // Matrices are stored in CSC format.
    const auto nz_rows = readDataset<uint32_t>(file, kNonZeroRows);
    const auto nz_data = readDataset<uint32_t>(file, kNonZeroData);
    const auto col_inds = readDataset<uint32_t>(file, kColumnIndices);
    if (nz_rows.size() != nz_data.size()) {
        stop(
            "Datasets \"%s\" and \"%s\" specify a different number of "
            "non-zero entries (%zu vs. %zu).",
            kNonZeroRows, kNonZeroData, nz_rows.size(), nz_data.size());
    }

    const std::vector<uint64_t> dims = getDimensions(file);
    matrix.init(MatrixMetadata{.nrow = static_cast<int>(dims[0]),
                               .ncol = static_cast<int>(dims[1]),
                               .nval = nz_data.size()});

    int col = 0;
    const int num_cols = col_inds.size() - 1;
    for (int i = 0; i < nz_rows.size(); i++) {
        while (col < num_cols && (col_inds[col + 1] <= i)) {
            col++;
        }
        if (col == num_cols) {
            break;
        }
        // Shift row and column by +1 for 1-based indexing.
        const MatrixData entry =
            MatrixData{.row = static_cast<int>(nz_rows[i]) + 1,
                       .col = col + 1,
                       .val = static_cast<int>(nz_data[i])};
        checkValid(matrix, entry);
        matrix.addEntry(entry);
    }
}

}  // namespace smallcount
