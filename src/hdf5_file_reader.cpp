#include "hdf5_file_reader.h"

#include <cstdint>
#include <string>
#include <vector>

#include "Rcpp.h"
#include "hdf5.h"
#include "sparse_matrix.h"
#include "tenx_file_params.h"

using namespace Rcpp;

namespace smallcount {
namespace {

// Name of the HDF5 group containing the matrix datasets.
std::string h5GroupName(const TenxFileParams &params) {
    return params.genome.has_value() ? *params.genome : "matrix";
}

// Name of the HDF5 dataset storing the non-zero row indices.
std::string indicesDataset(const TenxFileParams &params) {
    return h5GroupName(params) + "/indices";
}

// Name of the HDF5 dataset storing the non-zero data entries.
std::string dataDataset(const TenxFileParams &params) {
    return h5GroupName(params) + "/data";
}

// Name of the HDF5 dataset storing the non-zero column indices.
std::string indptrDataset(const TenxFileParams &params) {
    return h5GroupName(params) + "/indptr";
}

// Name of the HDF5 dataset storing the matrix dimensions.
std::string shapeDataset(const TenxFileParams &params) {
    return h5GroupName(params) + "/shape";
}

// Name of the HDF5 dataset storing the row names.
std::string featuresDataset(const TenxFileParams &params) {
    if (params.genome.has_value()) {
        std::string dataset =
            params.use_id_row_names ? "/genes" : "/gene_names";
        return h5GroupName(params) + dataset;
    }
    std::string dataset = params.use_id_row_names ? "/id" : "/name";
    return h5GroupName(params) + "/features" + dataset;
}

// Name of the HDF5 dataset storing the column names.
std::string barcodesDataset(const TenxFileParams &params) {
    return h5GroupName(params) + "/barcodes";
}

// Helper function to read HDF5 datasets of different types.
template <typename T>
std::vector<T> readData(hid_t dataset, hsize_t num_entries);
template <>
std::vector<uint32_t> readData<uint32_t>(hid_t dataset, hsize_t num_entries) {
    std::vector<uint32_t> data(num_entries);
    H5Dread(dataset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            data.data());
    return data;
}
template <>
std::vector<uint64_t> readData<uint64_t>(hid_t dataset, hsize_t num_entries) {
    std::vector<uint64_t> data(num_entries);
    H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            data.data());
    return data;
}
template <>
std::vector<std::string> readData<std::string>(hid_t dataset,
                                               hsize_t num_entries) {
    hid_t data_type = H5Dget_type(dataset);
    size_t str_size = H5Tget_size(data_type) + 1;

    hid_t mem_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(mem_type, str_size);

    char *buffer = new char[num_entries * str_size];
    H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    std::vector<std::string> data(num_entries);
    for (size_t i = 0; i < num_entries; i++) {
        std::string str(buffer + i * str_size, str_size);
        str.erase(str.find('\0'));
        data[i] = str;
    }

    delete[] buffer;
    H5Tclose(mem_type);
    H5Tclose(data_type);
    return data;
}

// Reads the contents of an HDF5 dataset into a vector of type T.
template <typename T>
std::vector<T> readDataset(hid_t file, const std::string &dataset_name) {
    hid_t dataset = H5Dopen2(file, dataset_name.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);

    hsize_t num_entries;
    H5Sget_simple_extent_dims(dataspace, &num_entries, NULL);
    std::vector<T> data = readData<T>(dataset, num_entries);

    H5Sclose(dataspace);
    H5Dclose(dataset);
    return data;
}

}  // namespace

SvtSparseMatrix Hdf5FileReader::read(hid_t file, const TenxFileParams &params) {
    // Check if HDF5 group exists.
    const std::string genome = h5GroupName(params);
    if (H5Lexists(file, genome.c_str(), H5P_DEFAULT) <= 0) {
        stop("Group '%s' not found in HDF5 file", genome.c_str());
    }

    // Read non-zero matrix entries in CSC format.
    const std::string indices_dataset = indicesDataset(params);
    const std::string data_dataset = dataDataset(params);
    const std::string indptr_dataset = indptrDataset(params);

    // Check that datasets exist.
    if (H5Lexists(file, indices_dataset.c_str(), H5P_DEFAULT) <= 0) {
        stop("Dataset '%s' not found in HDF5 file", indices_dataset.c_str());
    }
    if (H5Lexists(file, data_dataset.c_str(), H5P_DEFAULT) <= 0) {
        stop("Dataset '%s' not found in HDF5 file", data_dataset.c_str());
    }
    if (H5Lexists(file, indptr_dataset.c_str(), H5P_DEFAULT) <= 0) {
        stop("Dataset '%s' not found in HDF5 file", indptr_dataset.c_str());
    }

    const auto nz_rows = readDataset<uint32_t>(file, indices_dataset);
    const auto nz_data = readDataset<uint32_t>(file, data_dataset);
    const auto col_inds = readDataset<uint32_t>(file, indptr_dataset);
    if (nz_rows.size() != nz_data.size()) {
        stop(
            "Inconsistent HDF5 dataset sizes. Datasets \"%s\" and \"%s\" "
            "specify a different number of non-zero entries (%zu vs. %zu).",
            indices_dataset, data_dataset, nz_rows.size(), nz_data.size());
    }

    // Read matrix dimensions.
    const std::string shape_dataset = shapeDataset(params);
    if (H5Lexists(file, shape_dataset.c_str(), H5P_DEFAULT) <= 0) {
        stop("Dataset '%s' not found in HDF5 file", shape_dataset.c_str());
    }
    const auto dims = readDataset<uint64_t>(file, shape_dataset);
    if (dims.size() != 2) {
        stop(
            "Invalid matrix dimensions. Dataset \"%s\" has %zu entries "
            "(expected 2).",
            shape_dataset, dims.size());
    }

    // Read row and column names.
    const std::string features_dataset = featuresDataset(params);
    if (H5Lexists(file, features_dataset.c_str(), H5P_DEFAULT) <= 0) {
        stop("Dataset '%s' not found in HDF5 file", features_dataset.c_str());
    }
    auto row_names = readDataset<std::string>(file, features_dataset);
    if (row_names.size() != dims[0]) {
        warning(
            "Datasets \"%s\" and \"%s\" specify a different number of rows "
            "(%zu vs. %zu).",
            shape_dataset, features_dataset, dims[0], row_names.size());
    }
    std::vector<std::string> col_names{};
    if (params.use_barcode_col_names) {
        const std::string barcodes_dataset = barcodesDataset(params);
        if (H5Lexists(file, barcodes_dataset.c_str(), H5P_DEFAULT) <= 0) {
            stop("Dataset '%s' not found in HDF5 file",
                 barcodes_dataset.c_str());
        }
        col_names = readDataset<std::string>(file, barcodes_dataset);
        if (col_names.size() != dims[1]) {
            warning(
                "Datasets \"%s\" and \"%s\" specify a different number of "
                "columns (%zu vs. %zu).",
                shape_dataset, barcodes_dataset, dims[1], col_names.size());
        }
    }

    // Initialize SVT matrix with dims[1] columns, each containing 2 vectors
    // with row and value information.
    Svt svt(dims[1], SvtEntry(2));

    int col = 0;
    const int num_cols = col_inds.size() - 1;
    for (int i = 0; i < nz_rows.size(); i++) {
        while (col < num_cols && (col_inds[col + 1] <= i)) {
            col++;
        }
        if (col == num_cols) {
            break;
        }
        svt[col][kSvtRowInd].emplace_back(static_cast<int>(nz_rows[i]));
        svt[col][kSvtValInd].emplace_back(static_cast<int>(nz_data[i]));
    }

    return SvtSparseMatrix(std::move(svt),
                           MatrixMetadata{.nrow = static_cast<int>(dims[0]),
                                          .ncol = static_cast<int>(dims[1]),
                                          .nval = nz_data.size(),
                                          .row_names = std::move(row_names),
                                          .col_names = std::move(col_names)});
}

}  // namespace smallcount
