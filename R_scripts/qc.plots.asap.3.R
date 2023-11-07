##################################################
## Project: ASAP
## Script purpose: QC plots v3
## Date: 2023 October 18
## Author: Vincent Gardeux (vincent.gardeux@epfl.ch)
##################################################

# Parameters handling
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# Libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2)) # for theme
suppressPackageStartupMessages(library(plotly)) # for ggplotly
suppressPackageStartupMessages(source("hdf5_lib.R"))

# Functions
serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty=TRUE, digits = 3)
}

# Arguments
input_loom <- args[1]
input_dataset_path <- args[2]
output_dir <- args[3]

#input_loom <- "grrpvn_parsing_output.loom"
#input_dataset_path <- "/matrix"
#output_dir <- "./"

# Parameters
set.seed(42)
data.warnings <- NULL
time_idle <- 0
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  if(!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")
}

# Error case: Loom file does not exist
if(!file.exists(input_loom)) error.json(paste0("This file: '", input_loom, "', does not exist!"))

## Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
data.loom <- open_with_lock(input_loom, "r")
data.matrix <- fetch_dataset(data.loom, input_dataset_path, transpose = T) # If run on dimension reduction (like PCA), then do not transpose the matrix
data.mt.content <- fetch_dataset(data.loom, "/col_attrs/_Mitochondrial_Content", transpose = F)
close_all()

# Error case: Path in Loom file does not exist
if(is.null(data.matrix)) error.json(paste0("This file: '", input_loom, "', does not contain any dataset at path '", input_dataset_path, "'!"))

## Create Seurat object
data.seurat <- CreateSeuratObject(counts = data.matrix, min.cells = 0, min.features = 0) # Create our Seurat object using our data matrix (no filtering)

# Add metadata
data.seurat@meta.data[["percent.mt"]] <- data.mt.content

## Save some RAM
rm(data.matrix)
rm(data.loom)
rm(data.mt.content)

## Generate plots
p <- VlnPlot(data.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) & theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(plot = p, bg = "white", filename = paste0(output_dir,"qc.1.seurat.png"), width = 10, height = 5)
ggsave(plot = p, bg = "white", filename = paste0(output_dir,"qc.1.seurat.pdf"), width = 10, height = 5)
write(serialize(ggplotly(p)), file = paste0(output_dir,"qc.1.seurat.json"), append=F)

plot1 <- FeatureScatter(data.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2
ggsave(plot = p, bg = "white", filename = paste0(output_dir,"qc.2.seurat.png"), width = 10, height = 5)
ggsave(plot = p, bg = "white", filename = paste0(output_dir,"qc.2.seurat.pdf"), width = 10, height = 5)
write(serialize(ggplotly(p)), file = paste0(output_dir,"qc.2.seurat.json"), append=F)

# Generate default JSON file
stats <- list()
stats$time_idle = time_idle
stats$nber_rows = nrow(data.seurat)
stats$nber_cols = ncol(data.seurat)
if(!is.null(data.warnings)) stats$warnings = data.warnings
stats$plots <- list(paste0(output_dir,"qc.1.seurat.png"), paste0(output_dir,"qc.1.seurat.pdf"), paste0(output_dir,"qc.1.seurat.json"), paste0(output_dir,"qc.2.seurat.png"), paste0(output_dir,"qc.2.seurat.pdf"), paste0(output_dir,"qc.2.seurat.json"))
if(exists('output_dir') & !is.null(output_dir) & !is.na(output_dir)){
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"))
} else {
  cat(toJSON(stats, method="C", auto_unbox=T, digits = NA))
}

