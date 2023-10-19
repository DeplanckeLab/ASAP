##################################################
## Project: ASAP
## Script purpose: Cell Filtering v3
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
output_loom <- args[3]
param_nCount_RNA <- args[4] # UMI / reads, default = 0, in [0, +Inf]
if(is.null(param_nCount_RNA) || is.na(param_nCount_RNA) || param_nCount_RNA == "") param_nCount_RNA <- 0
param_nFeature_RNA <- as.numeric(args[5]) # detected genes, default = 0, in [0, +Inf]
if(is.null(param_nFeature_RNA) || is.na(param_nFeature_RNA) || param_nFeature_RNA == "") param_nFeature_RNA <- 0
param_percent.prot <- as.numeric(args[6]) # Percent protein coding genes, default = 0, in [0, 100]
if(is.null(param_percent.prot) || is.na(param_percent.prot) || param_percent.prot == "") param_percent.prot <- 0
param_percent.ribo <- as.numeric(args[7]) # Percent ribosomal genes, default = 100, in [0, 100]
if(is.null(param_percent.ribo) || is.na(param_percent.ribo) || param_percent.ribo == "") param_percent.ribo <- 100
param_percent.mt <- as.numeric(args[8]) # Percent mitochondrial genes, default = 100, in [0, 100]
if(is.null(param_percent.mt) || is.na(param_percent.mt) || param_percent.mt == "") param_percent.mt <- 100
data.warnings <- NULL
time_idle <- 0

#input_loom <- "grrpvn_parsing_output.loom"
#input_dataset_path <- "/matrix"
#output_loom <- "grrpvn_cell_filtering_output.loom"
#param_nCount_RNA <- 0 # UMI / reads, default = 0, in [0, +Inf]
#param_nFeature_RNA <- 0 # detected genes, default = 0, in [0, +Inf]
#param_percent.prot <- 95 # Percent protein coding genes, default = 0, in [0, 100]
#param_percent.ribo <- 100 # Percent ribosomal genes, default = 100, in [0, 100]
#param_percent.mt <- 100 # Percent mitochondrial genes, default = 100, in [0, 100]

# Error case: Loom file does not exist
if(!file.exists(input_loom)) error.json(paste0("This file: '", input_loom, "', does not exist!"))

## Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
list_row_attrs <- list()
list_col_attrs <- list()
data.loom <- open_with_lock(input_loom, "r")
data.matrix <- fetch_dataset(data.loom, input_dataset_path, transpose = T) # If run on dimension reduction (like PCA), then do not transpose the matrix
for(rowattr in names(data.loom$row_attrs)){
  list_row_attrs[[rowattr]] <- data.loom$row_attrs[[rowattr]]
}
for(colattr in names(data.loom$col_attrs)){
  list_col_attrs[[colattr]] <- data.loom$col_attrs[[colattr]]
}
# Should not have any layers or graph with the new architecture => SKIPPED
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.matrix)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", input_dataset_path, "'!"))

## Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.matrix, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

## Adding all cell metadata (we cannot add the gene metadata in the Seurat object?)
for(colattr in names(list_col_attrs)) {
  data.seurat@meta.data[[colattr]] <- list_col_attrs[[colattr]]
}
data.seurat@meta.data[["percent.mt"]] <- list_col_attrs[["_Mitochondrial_Content"]] # Should exist
data.seurat@meta.data[["percent.ribo"]] <- list_col_attrs[["_Ribosomal_Content"]] # Should exist
data.seurat@meta.data[["percent.prot"]] <- list_col_attrs[["_Protein_Coding_Content"]] # Should exist

## Save some RAM
rm(data.matrix)
rm(data.loom)

### Filter the matrix
data.seurat <- subset(data.seurat, subset = nCount_RNA >= param_nCount_RNA & nFeature_RNA >= param_nFeature_RNA & percent.prot >= param_percent.prot & percent.mt <= param_percent.mt & percent.ribo <= param_percent.ribo, return.null = T)
if(is.null(data.seurat)) error.json("Filtering failed. No more cells.")

## Filter the non-expressed genes
to_filter_genes <- (rowSums(data.seurat) > 0)
data.seurat <- data.seurat[to_filter_genes, ]

## Filter corresponding metadata
for(rowattr in names(list_row_attrs)){
  list_row_attrs[[rowattr]] <- list_row_attrs[[rowattr]][to_filter_genes]
}

# Prepare the JSON
stats <- list()
stats$metadata = list()
index <- 1

# Writing the dataset in the loom file
data.loom <- create_with_lock(output_loom)
# Add main matrix
add_matrix_dataset(handle = data.loom, dataset_path = input_dataset_path, dataset_object = t(as.matrix(data.seurat@assays$RNA@counts)))
stats$metadata[[index]] <- list(name = input_dataset_path, on = "EXPRESSION_MATRIX", type = "NUMERIC", nber_rows = nrow(data.seurat), nber_cols = ncol(data.seurat), dataset_size = get_dataset_size(data.loom, input_dataset_path))
index <- index + 1

# Add default groups
h5createGroup(file = data.loom, group = "/attrs")
h5createGroup(file = data.loom, group = "/col_attrs")
h5createGroup(file = data.loom, group = "/col_graphs")
h5createGroup(file = data.loom, group = "/layers")
h5createGroup(file = data.loom, group = "/row_attrs")
h5createGroup(file = data.loom, group = "/row_graphs")

# Add LOOM version
add_single_value(handle = data.loom, value_path = "/attrs/LOOM_SPEC_VERSION", value = "3.0.0", storage.mode_param = "character")

# Add row metadata
for(rowattr in names(list_row_attrs)){
  storage.mode <- typeof(list_row_attrs[[rowattr]])
  type.mode <- "STRING"
  if(storage.mode %in% c("integer", "double", "float", "long")) type.mode <- "NUMERIC"
  if(rowattr %in% c("_Biotypes",  "_Chromosomes", "variable_features")) type.mode <- "CATEGORICAL"
  # TODO: Handle properly the CATEGORICAL vs TEXT values
  add_array_dataset(handle = data.loom, dataset_path = paste0("/row_attrs/", rowattr), storage.mode_param = storage.mode, dataset_object = list_row_attrs[[rowattr]])
  stats$metadata[[index]] = list(name = paste0("/row_attrs/", rowattr), on = "GENE", type = type.mode, nber_rows = nrow(data.seurat), nber_cols = 1, dataset_size = get_dataset_size(data.loom, paste0("/row_attrs/", rowattr)))
  index <- index + 1
}

# Add col metadata
for(colattr in names(list_col_attrs)) {
  # TODO: Here we should probably recompute some data. Such as _Depth. Instead of just filtering it.
  storage.mode <- typeof(data.seurat@meta.data[[colattr]])
  type.mode <- "STRING"
  if(storage.mode %in% c("integer", "double", "float", "long")) type.mode <- "NUMERIC"
  add_array_dataset(handle = data.loom, dataset_path = paste0("/col_attrs/", colattr), storage.mode_param = storage.mode, dataset_object = data.seurat@meta.data[[colattr]])
  stats$metadata[[index]] = list(name = paste0("/col_attrs/", colattr), on = "CELL", type = type.mode, nber_rows = 1, nber_cols = ncol(data.seurat), dataset_size = get_dataset_size(data.loom, paste0("/col_attrs/", colattr)))
  index <- index + 1
}

close_all()

# Generate default JSON file
stats$time_idle = time_idle
stats$nber_rows = nrow(data.seurat)
stats$nber_cols = ncol(data.seurat)
if(!is.null(data.warnings)) stats$warnings = data.warnings
cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))

