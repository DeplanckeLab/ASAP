##################################################
## Project: ASAP
## Script purpose: UMAP calculation v3
## Date: 2023 October 19
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

### Parameters handling
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

## Libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(source("hdf5_lib.R"))

## Functions
serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty=TRUE, digits = 3)
}

# Arguments
input_loom <- args[1]
raw_dataset_path <- args[2]
pca_dataset_path <- args[3]
output_dataset_path <- args[4]
output_dir <- args[5]
n_pcs <- as.numeric(args[6]) # Number of PCs of PCA to use for computing UMAP

#input_loom <- "grrpvn_parsing_output.loom"
#raw_dataset_path <- "/matrix"
#pca_dataset_path <- "/col_attrs/PCA_seurat"
#output_dataset_path <- "/col_attrs/UMAP_seurat"
#n_pcs <- 50 # Default

# Parameters
set.seed(42)
data.warnings <- NULL
time_idle <- 0
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  if(!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")
}

# Error case: Loom file does not exist
if(!file.exists(input_loom)) error.json(paste0("This file: '", input_loom, "', does not exist!"))

# Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
data.loom <- open_with_lock(input_loom, "r")
data.raw <- fetch_dataset(data.loom, raw_dataset_path, transpose = T)
data.pca <- fetch_matrix(data.loom, pca_dataset_path, transpose = T)
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.raw)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", raw_dataset_path, "'!"))
if(is.null(data.pca)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", pca_dataset_path, "'!"))

# Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.raw, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

# Import PCA
colnames(data.pca) <- paste0('PC_', 1:ncol(data.pca))
rownames(data.pca) <- colnames(data.seurat)
data.seurat[['pca']] <- CreateDimReducObject(embeddings = data.pca, key = 'PC_', assay = 'RNA')

# Save some RAM
rm(data.raw)
rm(data.loom)
rm(data.pca)

# Run UMAP
suppressWarnings(data.seurat <- RunUMAP(data.seurat, reduction = "pca", dims = 1:n_pcs, umap.method = "uwot", verbose = F))

# Writing the UMAP in the loom file
data.loom <- open_with_lock(input_loom, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_dataset_path, storage.mode_param = "double", dataset_object = t(data.seurat@reductions$umap@cell.embeddings))
datasetSize <- get_dataset_size(data.loom, output_dataset_path)
close_all()

# Generate default JSON file
stats <- list()
stats$time_idle = time_idle
if(!is.null(data.warnings)) stats$warnings = data.warnings
stats$metadata = list(list(name = output_dataset_path, on = "CELL", type = "NUMERIC", nber_rows = 2, nber_cols = ncol(data.seurat), dataset_size = datasetSize))
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"))
} else {
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))
}
