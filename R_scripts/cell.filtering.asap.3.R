##################################################
## Project: ASAP
## Script purpose: Cell Filtering v3
## Date: 2023 October 12
## Updated: 2023 October 20
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Parameters handling
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# Libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(source("hdf5_lib.R"))

# Arguments
input_loom <- args[1]
input_dataset_path <- args[2]
output_loom <- args[3]
output_dir <- args[4]
param_nCount_RNA <- args[5] # UMI / reads, default = 0, in [0, +Inf]
if(is.null(param_nCount_RNA) || is.na(param_nCount_RNA) || param_nCount_RNA == "") param_nCount_RNA <- 0
param_nFeature_RNA <- as.numeric(args[6]) # detected genes, default = 0, in [0, +Inf]
if(is.null(param_nFeature_RNA) || is.na(param_nFeature_RNA) || param_nFeature_RNA == "") param_nFeature_RNA <- 0
param_percent.prot <- as.numeric(args[7]) # Percent protein coding genes, default = 0, in [0, 100]
if(is.null(param_percent.prot) || is.na(param_percent.prot) || param_percent.prot == "") param_percent.prot <- 0
param_percent.ribo <- as.numeric(args[8]) # Percent ribosomal genes, default = 100, in [0, 100]
if(is.null(param_percent.ribo) || is.na(param_percent.ribo) || param_percent.ribo == "") param_percent.ribo <- 100
param_percent.mt <- as.numeric(args[9]) # Percent mitochondrial genes, default = 100, in [0, 100]
if(is.null(param_percent.mt) || is.na(param_percent.mt) || param_percent.mt == "") param_percent.mt <- 100

# Parameters
set.seed(42)
data.warnings <- NULL
time_idle <- 0
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  if(!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")
}

#input_loom <- "/data/gardeux/f9w7cz_parsing_output.loom" #grrpvn_parsing_output.loom
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
list_layers <- list()
data.loom <- open_with_lock(input_loom, "r")
data.matrix <- fetch_dataset(data.loom, input_dataset_path, transpose = T) # If run on dimension reduction (like PCA), then do not transpose the matrix
for(rowattr in names(data.loom$row_attrs)){
  list_row_attrs[[rowattr]] <- data.loom$row_attrs[[rowattr]]
}
for(colattr in names(data.loom$col_attrs)){
  list_col_attrs[[colattr]] <- data.loom$col_attrs[[colattr]]
}
for(layer in names(data.loom$layers)){
  list_layers[[layer]] <- data.loom$layers[[layer]]
}
# For now, we don't care about the graphs, or attrs
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.matrix)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", input_dataset_path, "'!"))

## Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.matrix, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

## Adding main cell metadata (we cannot add the gene metadata in the Seurat object?)
data.seurat@meta.data[["percent.mt"]] <- list_col_attrs[["_Mitochondrial_Content"]] # Should exist
data.seurat@meta.data[["percent.ribo"]] <- list_col_attrs[["_Ribosomal_Content"]] # Should exist
data.seurat@meta.data[["percent.prot"]] <- list_col_attrs[["_Protein_Coding_Content"]] # Should exist

## Save some RAM
rm(data.matrix)
rm(data.loom)

### Filter the matrix
before_filtering <- colnames(data.seurat)
data.seurat <- subset(data.seurat, subset = nCount_RNA >= param_nCount_RNA & nFeature_RNA >= param_nFeature_RNA & percent.prot >= param_percent.prot & percent.mt <= param_percent.mt & percent.ribo <= param_percent.ribo, return.null = T)
if(is.null(data.seurat)) error.json("Filtering failed. No more cells.")
after_filtering <- colnames(data.seurat)

## Filter the non-expressed genes
to_filter_genes <- (rowSums(data.seurat) > 0)
data.seurat <- data.seurat[to_filter_genes, ]

## Filter row metadata
for(rowattr in names(list_row_attrs)){
  if(length(dim(list_row_attrs[[rowattr]])) == 1){
    list_row_attrs[[rowattr]] <- list_row_attrs[[rowattr]][to_filter_genes]
  } else if(length(dim(list_row_attrs[[rowattr]])) == 2){
    # DE results for example?
    list_row_attrs[[rowattr]] <- list_row_attrs[[rowattr]][, to_filter_genes]
  } else {
    # What is that?
  }
}

## Filter col metadata
to_filter_cells <- before_filtering %in% after_filtering 
for(colattr in names(list_col_attrs)){
  if(length(dim(list_col_attrs[[colattr]])) == 1){
    list_col_attrs[[colattr]] <- list_col_attrs[[colattr]][to_filter_cells]
  } else if(length(dim(list_col_attrs[[colattr]])) == 2){
    # Embeddings for example?
    list_col_attrs[[colattr]] <- list_col_attrs[[colattr]][, to_filter_cells]
  } else {
    # What is that?
  }
}

## Filter layers
for(layer in names(list_layers)){
  if(length(dim(list_layers[[layer]])) == 2){
    list_layers[[layer]] <- list_layers[[layer]][to_filter_cells, to_filter_genes]
  } else {
    # What is that?
  }
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
  # Check type
  storage.mode <- typeof(list_row_attrs[[rowattr]])
  type.mode <- "STRING"
  if(storage.mode %in% c("integer", "double", "float", "long")) type.mode <- "NUMERIC"
  if(rowattr %in% c("_Biotypes",  "_Chromosomes", "variable_features")) type.mode <- "CATEGORICAL"
  # TODO: Handle properly the CATEGORICAL vs TEXT values
  
  # Check dimension
  ncol <- 1
  if(length(dim(list_row_attrs[[rowattr]])) == 2) ncol <- dim(list_row_attrs[[rowattr]])[1]
  
  # Add to Loom
  if(ncol == 1) {
    add_array_dataset(handle = data.loom, dataset_path = paste0("/row_attrs/", rowattr), storage.mode_param = storage.mode, dataset_object = list_row_attrs[[rowattr]])
  } else {
    add_matrix_dataset(handle = data.loom, dataset_path = paste0("/row_attrs/", rowattr), storage.mode_param = storage.mode, dataset_object = list_row_attrs[[rowattr]])
  }

  # Prepare JSON
  stats$metadata[[index]] = list(name = paste0("/row_attrs/", rowattr), on = "GENE", type = type.mode, nber_rows = nrow(data.seurat), nber_cols = ncol, dataset_size = get_dataset_size(data.loom, paste0("/row_attrs/", rowattr)))
  index <- index + 1
}

# Add col metadata
for(colattr in names(list_col_attrs)) {
  # Check type
  storage.mode <- typeof(list_col_attrs[[colattr]])
  type.mode <- "STRING"
  if(storage.mode %in% c("integer", "double", "float", "long")) type.mode <- "NUMERIC"
  # TODO: Handle properly the CATEGORICAL vs TEXT values
  
  # Check dimension
  nrow <- 1
  if(length(dim(list_col_attrs[[colattr]])) == 2) nrow <- dim(list_col_attrs[[colattr]])[1]
  
  # Add to Loom
  if(nrow == 1) {
    add_array_dataset(handle = data.loom, dataset_path = paste0("/col_attrs/", colattr), storage.mode_param = storage.mode, dataset_object = list_col_attrs[[colattr]])
  } else {
    add_matrix_dataset(handle = data.loom, dataset_path = paste0("/col_attrs/", colattr), storage.mode_param = storage.mode, dataset_object = list_col_attrs[[colattr]])
  }
  
  # Prepare JSON
  stats$metadata[[index]] = list(name = paste0("/col_attrs/", colattr), on = "CELL", type = type.mode, nber_rows = nrow, nber_cols = ncol(data.seurat), dataset_size = get_dataset_size(data.loom, paste0("/col_attrs/", colattr)))
  index <- index + 1
}

# Add layers
for(layer in names(list_layers)) {
  # Check type
  storage.mode <- typeof(list_layers[[layer]])
  type.mode <- "STRING"
  if(storage.mode %in% c("integer", "double", "float", "long")) type.mode <- "NUMERIC"
  # TODO: Handle properly the CATEGORICAL vs TEXT values
  
  # Check dimension
  if(length(dim(list_layers[[layer]])) == 2){
    nrow <- dim(list_layers[[layer]])[2]
    ncol <- dim(list_layers[[layer]])[1]
    
    # Add to Loom
    add_matrix_dataset(handle = data.loom, dataset_path = paste0("/layers/", layer), storage.mode_param = storage.mode, dataset_object = list_layers[[layer]])

    # Prepare JSON
    stats$metadata[[index]] = list(name = paste0("/layers/", colattr), on = "EXPRESSION_MATRIX", type = type.mode, nber_rows = nrow, nber_cols = ncol, dataset_size = get_dataset_size(data.loom, paste0("/layers/", layer)))
    index <- index + 1
  } else {
    # What is that?
  }
}

# Recompute standard stuff
add_array_dataset(handle = data.loom, dataset_path = "/row_attrs/_Sum", storage.mode_param = "double", dataset_object = rowSums(data.seurat))
add_array_dataset(handle = data.loom, dataset_path = "/col_attrs/_Depth", storage.mode_param = "double", dataset_object = colSums(data.seurat))
add_array_dataset(handle = data.loom, dataset_path = "/row_attrs/_Detected_Genes", storage.mode_param = "double", dataset_object = colSums(data.seurat@assays$RNA@counts > 0))

close_all()

# Generate default JSON file
stats$time_idle = time_idle
stats$nber_rows = nrow(data.seurat)
stats$nber_cols = ncol(data.seurat)
if(!is.null(data.warnings)) stats$warnings = data.warnings
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"))
} else {
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))
}
