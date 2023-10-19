##################################################
## Project: ASAP
## Script purpose: Normalization v3
## Date: 2023 November 12
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

### Parameters handling
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

## Libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(source("hdf5_lib.R"))

### Default Parameters
set.seed(42)
input_loom <- args[1]
input_dataset_path <- args[2]
output_dataset_path <- args[3]
data.warnings <- NULL
time_idle <- 0

#input_loom <- "grrpvn_parsing_output.loom"
#input_dataset_path <- "/matrix"
#output_dataset_path <- "/layers/normalized"

# Error case: Loom file does not exist
if(!file.exists(input_loom)) error.json(paste0("This file: '", input_loom, "', does not exist!"))
 
## Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
data.loom <- open_with_lock(input_loom, "r")
data.matrix <- fetch_dataset(data.loom, input_dataset_path, transpose = T) # If run on dimension reduction (like PCA), then do not transpose the matrix
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.matrix)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", input_dataset_path, "'!"))

## Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.matrix, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

## Save some RAM
rm(data.matrix)
rm(data.loom)

### Normalize the raw matrix
data.seurat <- NormalizeData(data.seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)

# Writing the dataset in the loom file
data.loom <- open_with_lock(input_loom, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_dataset_path, dataset_object = t(as.matrix(data.seurat@assays$RNA@data)))
datasetSize <- get_dataset_size(data.loom, output_dataset_path)
close_all()

# Generate default JSON file
stats <- list()
stats$time_idle = time_idle
stats$nber_rows = nrow(data.seurat)
stats$nber_cols = ncol(data.seurat)
if(!is.null(data.warnings)) stats$warnings = data.warnings
stats$metadata = list(list(name = output_dataset_path, on = "EXPRESSION_MATRIX", type = "NUMERIC", nber_rows = nrow(data.seurat), nber_cols = ncol(data.seurat), dataset_size = datasetSize))
cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))

