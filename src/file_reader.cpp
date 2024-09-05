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

void readCsvFile(const std::string &filepath, SparseMatrix &matrix) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        stop("Could not open file: %s", filepath);
    }
    CsvFileReader::read(file, matrix);
    file.close();
}

void readMtxFile(const std::string &filepath, SparseMatrix &matrix) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        stop("Could not open file: %s", filepath);
    }
    MtxFileReader::read(file, matrix);
    file.close();
}

void readHdf5File(const std::string &filepath, SparseMatrix &matrix) {
    hid_t file = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        stop("Could not open file: %s", filepath);
    }
    Hdf5FileReader::read(file, matrix);
    H5Fclose(file);
}

void validateArguments(const std::string &file_extension,
                       const std::string &rep, SparseMatrix *matrix) {
    if (supportedExtensions.find(file_extension) == supportedExtensions.end()) {
        stop("File extension %s not supported. Try *.csv, *.mtx, or *.h5.",
             file_extension);
    }
    if (matrix == nullptr) {
        stop(
            "Invalid matrix representation specified: \"%s\". Only \"coo\" "
            "and \"svt\" supported.",
            rep);
    }
}

}  // namespace

SEXP SparseMatrixFileReader::read(const std::string &filepath,
                                  const std::string &rep) {
    const std::string file_extension = get_extension(filepath);
    std::unique_ptr<SparseMatrix> matrix = SparseMatrix::create(rep);
    validateArguments(file_extension, rep, matrix.get());
    if (file_extension == kCsv) {
        readCsvFile(filepath, *matrix);
    } else if (file_extension == kMtx) {
        readMtxFile(filepath, *matrix);
    } else {
        readHdf5File(filepath, *matrix);
    }
    return matrix->toRcpp();
}

}  // namespace smallcount
