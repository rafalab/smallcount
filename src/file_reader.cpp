#include "file_reader.h"

#include <fstream>
#include <memory>
#include <string>
#include <unordered_set>

#include "Rcpp.h"
#include "csv_file_reader.h"
#include "hdf5.h"
#include "hdf5_file_reader.h"
#include "mtx_file_reader.h"
#include "sparse_matrix.h"
#include "tenx_file_params.h"

using namespace Rcpp;

namespace smallcount {
namespace {

static constexpr char kCsv[] = "csv";
static constexpr char kMtx[] = "mtx";
static constexpr char kHdf5[] = "h5";

const std::unordered_set<std::string> supportedExtensions = {kCsv, kMtx, kHdf5};

inline std::string get_extension(const std::string &filepath) {
    return filepath.substr(filepath.find_last_of(".") + 1);
}

std::ifstream openFile(const std::string &filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        stop("Could not open file: %s", filepath);
    }
    return file;
}

SEXP readCsvFile(const std::string &filepath) {
    std::ifstream file = openFile(filepath);
    SvtSparseMatrix matrix = CsvFileReader::read(file);
    file.close();

    return matrix.toRcpp();
}

SEXP readMtxFile(const std::string &filedir, const TenxFileParams &params) {
    std::ifstream matrix_file = openFile(filedir + "matrix.mtx");
    std::ifstream barcodes_file = openFile(filedir + "barcodes.tsv");
    const std::string features_filename =
        params.use_features_tsv ? "features.tsv" : "genes.tsv";
    std::ifstream features_file = openFile(filedir + features_filename);
    SvtSparseMatrix matrix =
        MtxFileReader::read(matrix_file, barcodes_file, features_file, params);
    matrix_file.close();
    barcodes_file.close();
    features_file.close();

    return matrix.toRcpp();
}

SEXP readHdf5File(const std::string &filepath, const TenxFileParams &params) {
    hid_t file = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        stop("Could not open file: %s", filepath);
    }
    SvtSparseMatrix matrix = Hdf5FileReader::read(file, params);
    H5Fclose(file);

    return matrix.toRcpp();
}

}  // namespace

SEXP SparseMatrixFileReader::read(const std::string &filepath,
                                  const TenxFileParams &params) {
    const std::string file_extension = get_extension(filepath);
    if (file_extension == kCsv) {
        return readCsvFile(filepath);
    } else if (file_extension == kHdf5) {
        return readHdf5File(filepath, params);
    }
    return readMtxFile(filepath, params);
}

}  // namespace smallcount
