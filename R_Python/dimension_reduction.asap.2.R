### ASAP Dim. Reduction script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
require(jsonlite)
source("hdf5_lib.R")

### Default Parameters
set.seed(42)
input_matrix_filename <- args[1]
output_dir <- args[2]
std_method_name <- args[3]
input_matrix_dataset <- args[4]
output_matrix_dataset <- args[5]
nber_dims <- as.numeric(args[6])
data.warnings <- NULL
time_idle <- 0

# UMAP
RunUMAP <- function(data.use, max.dim = 2, n_neighbors = 30, min_dist = 0.3, metric = "correlation") {
  if (!py_module_available(module = 'umap')) error.json("Cannot find UMAP, please install Python package through pip (e.g. pip install umap-learn).")
  umap_import <- import("umap")
  umap <- umap_import$UMAP(random_state=as.integer(x = 42), n_neighbors = as.integer(x = n_neighbors), n_components = as.integer(x = max.dim), metric = metric, min_dist = min_dist)
  umap_output <- umap$fit_transform(as.matrix(x = data.use))
  return(umap_output)
}

### Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
to_transpose <- T
if(startsWith(input_matrix_dataset, "/col_attrs")) to_transpose <- F
data.loom <- open_with_lock(input_matrix_filename, "r")
data.parsed <- fetch_dataset(data.loom, input_matrix_dataset, transpose = to_transpose) # If run on dimension reduction (like PCA), then do not transpose the matrix
close_all()
if(is.null(data.parsed)) error.json("This dataset does not exist in the Loom file")

### Preprocess the data
nbNoVar = nrow(data.parsed)
data.parsed <- data.parsed[apply(data.parsed, 1, var, na.rm=TRUE) != 0,] # Remove lines without variation
nbNoVar = nbNoVar - nrow(data.parsed)
message("Removed ", nbNoVar, " rows without variation.")
if(ncol(data.parsed) < nber_dims | nrow(data.parsed) < nber_dims) error.json(paste0("Your matrix is too small to run the dimension reduction (less dimensions that the nb dims requested: ",nber_dims,")"))

### Remove duplicates
#data.parsed <- data.parsed[,!duplicated(t(data.parsed))]

### Run dimension reduction std_method_name
if (std_method_name == "pca"){ # default []
  # Packages
  require(irlba)
  
  # Parameters
  isCenter = args[7]
  if(!is.null(isCenter) & !is.na(isCenter)) isCenter = (isCenter == "true")
  else isCenter = T
  
  isScale = args[8]
  if(!is.null(isScale) & !is.na(isScale)) isScale = (isScale == "true")
  else isScale = F
  
  tryCatch({
    data.out <<- prcomp_irlba(t(data.parsed), center = isCenter, scale. = isScale, n = nber_dims)$x
  }, error = function(err) {
    if(grepl("convergence criterion below machine epsilon", err)) error.json("[irlba ERROR] Convergence criterion is below machine epsilon")
    error.json(err)
  })
} else if (std_method_name == "tsne"){
  # Packages
  require(Rtsne)
  perp.max <- (ncol(data.parsed) - 1) / 3
  
  # Parameters
  perp <- args[7]
  if(!is.null(perp) & !is.na(perp)) perp <- as.numeric(perp)
  else perp <- perp.max
  if(perp > perp.max){
    data.warnings <<- list(name=paste0("Perplexity is too high. Set to ", perp.max, "."), description=paste0("Perplexity should not be over nbCell/3."))
    perp <- perp.max
  }
  
  thet <- args[8]
  if(!is.null(thet) & !is.na(thet)) thet <- as.numeric(thet)
  else thet <- 0.5
  
  tryCatch({
    data.tsne <<- Rtsne(t(data.parsed), dims = nber_dims, theta = thet, check.duplicates = F, perplexity = perp)
  }, error = function(err) {
    if(grepl("Remove duplicates", err$message)) error.json("[Rtsne ERROR] Remove duplicates before running t-SNE")
    error.json(err$message)
  })
  
  data.out <- data.tsne$Y
} else if (std_method_name == "umap"){
  # Packages
  require(reticulate)
  
  # Parameters
  min_dist <- args[7]
  if(!is.null(min_dist) & !is.na(min_dist)) min_dist <- as.numeric(min_dist)
  else min_dist <- 0.2

  n_neighbors <- args[8]
  if(!is.null(n_neighbors) & !is.na(n_neighbors)) n_neighbors <- as.numeric(n_neighbors)
  else n_neighbors <- 30
  if(n_neighbors > ncol(data.parsed)){
    n_neighbors <- ncol(data.parsed) - 1
    data.warnings <<- list(name=paste0("n_neighbors is larger than the dataset size; truncating to ", n_neighbors, "."), description=paste0("Truncated to nbCells - 1."))
  }
  
  metric <- args[9]
  if(is.null(metric) | is.na(metric)) metric <- "correlation"
  
  # Run UMAP
  tryCatch({
    data.out <<- RunUMAP(data.use = t(data.parsed), min_dist = min_dist, n_neighbors = n_neighbors, max.dim = nber_dims, metric = metric)
  }, error = function(err) {
    if(grepl("ValueError: Input contains NaN", err)) error.json("[UMAP ERROR] Input contains NaN, infinity or a value too large for dtype(float32)")
    error.json(err)
  })
} else error.json("This dimension reduction method is not implemented.")

# Open Loom in writing mode for writing results
data.loom <- open_with_lock(input_matrix_filename, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_matrix_dataset, dataset_object = t(data.out))
close_all()

# Generate default JSON file
stats <- list()
stats$time_idle <- time_idle
# Prepare metadata report
stats$metadata = list(list(name = output_matrix_dataset, on = "CELL", type = "NUMERIC", nber_cols = nrow(data.out), nber_rows = ncol(data.out)))
if(!is.null(data.warnings)) stats$warnings = list(data.warnings)
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
