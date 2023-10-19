### ASAP H5ad script
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
suppressPackageStartupMessages(require(jsonlite))
suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(SeuratDisk))
suppressPackageStartupMessages(require(hdf5r))

### Default Parameters
set.seed(42)
input_h5ad_path <- args[1]
output_dir <- args[2]
input_h5ad_path <- "/data/vincent/s_fca_biohub_gut_10x.h5ad"

Convert(source = "/data/vincent/s_fca_biohub_gut_10x.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat.object <- LoadH5Seurat("/data/vincent/s_fca_biohub_gut_10x.h5seurat")


### Open the .h5 file in read/write and add datasets
handle <- H5Fopen(name = input_h5ad_path, flags = "H5F_ACC_RDONLY")

sparse_format <- h5readAttributes(file = handle, name = "/X")$h5sparse_format
sparse_shape <- h5readAttributes(file = handle, name = "/X")$h5sparse_shape

matrix = h5read(handle,"/X")
sparseMatrix(i = matrix$indices, j = ep, p = , repr = "C")
h5createGroup(file = handle,"mtx")
h5write(file = handle, obj = barcode.names, name = "/mtx/barcodes", level = 2)
h5write(file = handle, obj = feature.names, name = "/mtx/genes", level = 2)
h5write(file = handle, obj = feature.names, name = "/mtx/gene_names", level = 2)
h5write(file = handle, obj = dim(mat), name = "/mtx/shape", level = 2)
h5write(file = handle, obj = mat@x, name = "/mtx/data", level = 2)
h5write(file = handle, obj = mat@i, name = "/mtx/indices", level = 2)
h5write(file = handle, obj = mat@p, name = "/mtx/indptr", level = 2)
h5closeAll()
