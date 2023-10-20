##################################################
## Project: ASAP
## Script purpose: Scaling and covariate removal v3
## Date: 2023 November 16
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
output_dataset_path <- args[4]
output_dir <- args[5]
vars_to_regress_param <- args[6]

#input_loom <- "grrpvn_parsing_output.loom"
#raw_dataset_path <- "/matrix"
#norm_dataset_path <- "/layers/normalized"
#output_dataset_path <- "/layers/scaled"
#vars_to_regress_param <- NULL
#vars_to_regress_param <- ""
#vars_to_regress_param <- "/col_attrs/_Mitochondrial_Content"
#vars_to_regress_param <- "/col_attrs/_Mitochondrial_Content,/col_attrs/_Depth"

vars_to_regress <- c()
if(!is.null(vars_to_regress_param) && !is.na(vars_to_regress_param) && vars_to_regress_param != "") vars_to_regress <- strsplit(vars_to_regress_param, split = ",")[[1]]
vars_to_regress_name <- c()

# Regress check
for(var in vars_to_regress){
  if(!startsWith(var, "/col_attrs/")) {
    error.json(paste0("This var to regress: '", var, "', should start with /col_attrs/"))
  } else {
    vars_to_regress_name <- c(vars_to_regress_name, gsub(pattern = "/col_attrs/",replacement = "", x = var))
  }
}

# Parameters
set.seed(42)
data.warnings <- NULL
time_idle <- 0
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  if(!endsWith(output_dir, "/")) output_dir <- output_dir + "/"
}

# Error case: Loom file does not exist
if(!file.exists(input_loom)) error.json(paste0("This file: '", input_loom, "', does not exist!"))

# Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
data.loom <- open_with_lock(input_loom, "r")
data.raw <- fetch_dataset(data.loom, raw_dataset_path, transpose = T) # If run on dimension reduction (like PCA), then do not transpose the matrix
data.norm <- fetch_matrix(data.loom, norm_dataset_path, transpose = T) # If run on dimension reduction (like PCA), then do not transpose the matrix
data.vars_to_regress <- list()
for(var in vars_to_regress_name) data.vars_to_regress[[var]] <- fetch_dataset(data.loom, paste0("/col_attrs/", var), transpose = F)
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.raw)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", raw_dataset_path, "'!"))
if(is.null(data.norm)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", norm_dataset_path, "'!"))

# Error case: Weird size error?
if(!all(dim(data.norm) == dim(data.raw))) error.json("Raw and Normalized datasets are not of the same size?")

# Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.raw, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

# Import Normalization
data.seurat@assays$RNA@data <- as(data.norm, "dgCMatrix")
colnames(data.seurat@assays$RNA@data) <- colnames(data.seurat@assays$RNA@counts)
rownames(data.seurat@assays$RNA@data) <- rownames(data.seurat@assays$RNA@counts)

# Add metadata
corrected_names <- c()
for(var in vars_to_regress_name) {
  varcorrected <- var
  if(startsWith(var, prefix = "_")) varcorrected <- substr(var, 2, nchar(var))
  data.seurat@meta.data[[varcorrected]] <- data.vars_to_regress[[var]]
  corrected_names <- c(corrected_names, varcorrected)
}

# Save some RAM
rm(data.raw)
rm(data.norm)
rm(data.loom)
rm(data.vars_to_regress)

# Scale data
all.genes <- rownames(data.seurat) # Here I want to scale all genes. By default, only HVG genes are scaled (speed up computation)
data.seurat <- ScaleData(data.seurat, features = all.genes, vars.to.regress = corrected_names, verbose = T)

# Writing the scaled matrix in the loom file
data.loom <- open_with_lock(input_loom, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_dataset_path, dataset_object = t(as.matrix(data.seurat@assays$RNA@scale.data)))
datasetSize <- get_dataset_size(data.loom, output_dataset_path)
close_all()

# Generate default JSON file
stats <- list()
stats$time_idle = time_idle
stats$nber_rows = nrow(data.seurat)
stats$nber_cols = ncol(data.seurat)
if(!is.null(data.warnings)) stats$warnings = data.warnings
stats$metadata = list(list(name = output_dataset_path, on = "EXPRESSION_MATRIX", type = "NUMERIC", nber_rows = nrow(data.seurat), nber_cols = ncol(data.seurat), dataset_size = datasetSize))
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"))
} else {
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))
}
