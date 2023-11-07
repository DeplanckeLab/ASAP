##################################################
## Project: ASAP
## Script purpose: PCA calculation v3
## Date: 2023 October 18
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
norm_dataset_path <- args[3]
scaled_dataset_path <- args[4]
variable_features_path <- args[5]
output_dataset_path <- args[6]
output_dir <- args[7]
n_pcs <- as.numeric(args[8]) # Number of PCs to compute

#input_loom <- "grrpvn_parsing_output.loom"
#raw_dataset_path <- "/matrix"
#norm_dataset_path <- "/layers/normalized"
#scaled_dataset_path <- "/layers/scaled"
#variable_features_path <- "/row_attrs/variable_features"
#output_dataset_path <- "/col_attrs/PCA_seurat"
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
data.norm <- fetch_matrix(data.loom, norm_dataset_path, transpose = T)
data.scaled <- fetch_matrix(data.loom, scaled_dataset_path, transpose = T)
data.variable_features <- fetch_dataset(data.loom, variable_features_path, transpose = F)
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.raw)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", raw_dataset_path, "'!"))
if(is.null(data.norm)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", norm_dataset_path, "'!"))
if(is.null(data.scaled)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", scaled_dataset_path, "'!"))
if(is.null(data.variable_features)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", variable_features_path, "'!"))

# Error case: Weird size error?
if(!all(dim(data.norm) == dim(data.raw))) error.json("Raw and Normalized datasets are not of the same size?")
if(!all(dim(data.scaled) == dim(data.raw))) error.json("Raw and Scaled datasets are not of the same size?")

# Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.raw, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

# Import Normalization
data.seurat@assays$RNA@data <- as(data.norm, "dgCMatrix")
colnames(data.seurat@assays$RNA@data) <- colnames(data.seurat@assays$RNA@counts)
rownames(data.seurat@assays$RNA@data) <- rownames(data.seurat@assays$RNA@counts)

# Import Scaling
data.seurat@assays$RNA@scale.data <- data.scaled
colnames(data.seurat@assays$RNA@scale.data) <- colnames(data.seurat@assays$RNA@counts)
rownames(data.seurat@assays$RNA@scale.data) <- rownames(data.seurat@assays$RNA@counts)

# Save some RAM
rm(data.raw)
rm(data.norm)
rm(data.scaled)
rm(data.loom)

# Run PCA
data.seurat <- RunPCA(data.seurat, features = rownames(data.seurat)[data.variable_features], npcs = n_pcs, verbose = F)

# Writing the PCA in the loom file
data.loom <- open_with_lock(input_loom, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_dataset_path, storage.mode_param = "double", dataset_object = t(data.seurat@reductions$pca@cell.embeddings))
datasetSize <- get_dataset_size(data.loom, output_dataset_path)
close_all()

# Generate default JSON file
stats <- list()
stats$time_idle = time_idle
if(!is.null(data.warnings)) stats$warnings = data.warnings
stats$metadata = list(list(name = output_dataset_path, on = "CELL", type = "NUMERIC", nber_rows = n_pcs, nber_cols = ncol(data.seurat), dataset_size = datasetSize))
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"))
} else {
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))
}
