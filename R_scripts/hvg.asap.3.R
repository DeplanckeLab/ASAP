##################################################
## Project: ASAP
## Script purpose: HVG calculation v3
## Date: 2023 October 13
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
method <- args[6] # vst, dispersion, mean.var.plot
n_features <- as.numeric(args[7]) # Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'

#input_loom <- "grrpvn_parsing_output.loom"
#raw_dataset_path <- "/matrix"
#norm_dataset_path <- "/layers/normalized"
#output_dataset_path <- "/row_attrs/variable_features"
#output_dir <- "./"
#n_features <- 2000 # Default
#method <- "mean.var.plot"

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

# Save some RAM
rm(data.raw)
rm(data.norm)
rm(data.loom)

# Find variable genes
data.seurat <- FindVariableFeatures(data.seurat, selection.method = method, nfeatures = n_features, verbose = F) # Here I arbitrarily take 2000 top HVG (default)
variableFeatures <- VariableFeatures(data.seurat)

# Plot variable gene graphics
p <- VariableFeaturePlot(data.seurat, raster = F)
ggsave(plot = p, bg = "white", filename = paste0(output_dir,"hvg.seurat.png"), width = 8, height = 5)
ggsave(plot = p, bg = "white", filename = paste0(output_dir,"hvg.seurat.pdf"), width = 8, height = 5)
write(serialize(ggplotly(p)), file = paste0(output_dir,"hvg.seurat.json"), append=F)

# Writing the hvg in the loom file
data.loom <- open_with_lock(input_loom, "r+")
add_array_dataset(handle = data.loom, dataset_path = output_dataset_path, storage.mode_param = "logical", dataset_object = rownames(data.seurat) %in% variableFeatures)
datasetSize <- get_dataset_size(data.loom, output_dataset_path)
close_all()

# Generate default JSON file
stats <- list()
stats$time_idle = time_idle
if(!is.null(data.warnings)) stats$warnings = data.warnings
stats$plots <- list(paste0(output_dir,"hvg.seurat.png"), paste0(output_dir,"hvg.seurat.pdf"), paste0(output_dir,"hvg.seurat.json"))
stats$metadata = list(list(name = output_dataset_path, on = "GENE", type = "DISCRETE", nber_rows = nrow(data.seurat), nber_cols = 1, dataset_size = datasetSize))
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"))
} else {
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))
}
