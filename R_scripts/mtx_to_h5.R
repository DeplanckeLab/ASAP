### ASAP mtxToH5 script
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
#
### Libraries
suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(hdf5r))
#
### Default Parameters
set.seed(42)
input_folder <- args[1]
output_h5_path <- args[2]
#
### Generate correct paths
barcode.path <- paste0(input_folder, "/barcodes.tsv")
features.path <-  paste0(input_folder, "/features.tsv")
if (!file.exists(features.path)){
features.path <-  paste0(input_folder, "/genes.tsv")
}
matrix.path <-  paste0(input_folder, "/matrix.mtx")
#
### Read sparse Matrix (.mtx)
mat <- readMM(file = matrix.path)
#
### Transform to dgCMatrix format (hdf5 10x handled)
mat <- as(mat, "dgCMatrix")
#
### Read features and barcodes
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)$V1
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)$V1
#
### Remove output file if already there
file.remove(output_h5_path)
#
### Create .h5 file
h5createFile(output_h5_path)
#
### Open the .h5 file in read/write and add datasets
handle <- H5Fopen(name = output_h5_path, flags = "H5F_ACC_RDWR")
h5createGroup(file = handle,"mtx")
h5write(file = handle, obj = barcode.names, name = "/mtx/barcodes", level = 2)
h5write(file = handle, obj = feature.names, name = "/mtx/genes", level = 2)
h5write(file = handle, obj = feature.names, name = "/mtx/gene_names", level = 2)
h5write(file = handle, obj = dim(mat), name = "/mtx/shape", level = 2)
h5write(file = handle, obj = mat@x, name = "/mtx/data", level = 2)
h5write(file = handle, obj = mat@i, name = "/mtx/indices", level = 2)
h5write(file = handle, obj = mat@p, name = "/mtx/indptr", level = 2)
h5closeAll()