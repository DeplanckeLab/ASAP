### ASAP Normalization script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
require(jsonlite)
require(loomR) # For handling Loom files

### Default Parameters
set.seed(42)
input_matrix_filename <- args[1]
output_dir <- args[2]
std_method_name <- args[3]
output_matrix_dataset <- args[4]
to.transpose <- T # Do I need to transpose the matrix in the Loom
toKeep <- NULL
data.parsed <- data.frame()
data.gene.length <- NULL
data.sum <- NULL
data.ercc <- NULL
nber_genes <- 0
nber_cells <- 0
nber_zeros <- 0

### Loom Handling
chunk.dims_param = c(256, 256)

#### Functions ####
LogNorm <- function(data, scale_factor, display_progress = TRUE) {
  .Call('_Seurat_LogNorm', PACKAGE = 'Seurat', data, scale_factor, display_progress)
}

FastRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastRowScale', PACKAGE = 'Seurat', mat, scale, center, scale_max, display_progress)
}

# Because there is a bug with the name args and I don't need Regressing out at this stage
ScaleData.loom <- function(object, name = 'layers/scale_data', do.scale = TRUE, do.center = TRUE, scale.max = 10, chunk.size = 1000, normalized.data = 'layers/norm_data', overwrite = FALSE, display.progress = TRUE, dtype = h5types$float, chunk.dims = c(256, 256)) {
  # Run Scaling
  object$apply(
    name = name,
    FUN = function(mat) {
      return(t(x = FastRowScale(
        mat = t(x = mat),
        scale = do.scale,
        center = do.center,
        scale_max = scale.max,
        display_progress = F
      )))
    },
    MARGIN = 1,
    chunk.size = chunk.size,
    dataset.use = normalized.data,
    overwrite = overwrite,
    display.progress = display.progress,
    dtype = dtype,
    chunk.dims = chunk.dims
  )
  # Clean up
  gc(verbose = FALSE)
  invisible(x = object)
}

# Because there is a bug with the name args
NormalizeData.loom <- function(object, scale.factor = 1e4, chunk.size = 1000, name = 'layers/norm_data', dataset.use = 'matrix', display.progress = TRUE, overwrite = FALSE, dtype = h5types$float, chunk.dims = c(256, 256)) {
  # Run Normalization
  object$apply(
    name = name,
    FUN = function(mat) {
      return(t(x = as.matrix(x = LogNorm(
        data = Matrix(data = t(x = mat), sparse = TRUE),
        scale_factor = scale.factor,
        display_progress = FALSE
      ))))
    },
    MARGIN = 2,
    chunk.size = chunk.size,
    dataset.use = dataset.use,
    overwrite = overwrite,
    display.progress = display.progress,
    dtype = dtype,
    chunk.dims = chunk.dims
  )
  # Clean up
  gc(verbose = FALSE)
  invisible(x = object)
}

error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output_dir,"/output.json"), append=F)
  stop(displayed)
}

# Open the Loom file while handling potential locking
open_with_lock <- function(loom_filename, mode) {
  repeat{ # Handle the lock of the file
    isLocked <- F
    tryCatch({
      data.loom <- connect(filename = loom_filename, mode = mode)
    }, error = function(err) {
      if(grepl("unable to lock file", err$message)) isLocked <<- T
      else error.json(err$message)
    })
    if(!isLocked) return(data.loom)
    else {
      message("Sleeping 1sec for file lock....")
      time_idle <<- time_idle + 1
      Sys.sleep(1)
    }
  }
}

# Copy the old Loom (we need to lock the file, so nobody writes in it while we copy)
data.loom <- open_with_lock(input_matrix_filename, "r") # Block access to the file (nobody can write while we copy)
if(std_method_name != "seurat") { 
  data.parsed <- as.data.frame(t(data.loom$matrix[, ])) # t() because loomR returns the t() of the correct matrix we want
} else file.copy(from = input_matrix_filename, to = paste0(output_dir, "/output.loom"), overwrite = T, recursive = F, copy.mode = T, copy.date = F)
if(std_method_name %in% c("tpm", "rpkm")) {
  data.gene.length <- data.loom$row.attrs$`_SumExonLength`[]
  toKeep <- data.gene.length > 0
  data.gene.length[!toKeep] = NA
}
if(std_method_name %in% c("cpm", "tpm", "rpkm")) data.sum <- data.loom$col.attrs$`_Depth`[]
if(std_method_name == "ruvg") {
  if(!data.loom$exists("/col_attrs/_ERCCs")) error.json(displayed = "RUVg method can only be ran with ERCCs (spike-ins) used as 'negative control genes'")
  data.ercc <- as.data.frame(data.loom$col.attrs$`_ERCCs`[, ])
  if(nrow(data.ercc) == 0) error.json(displayed = "RUVg method can only be ran with ERCCs (spike-ins) used as 'negative control genes'")
}
data.loom$close_all()

### Functions
error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output_dir,"/output.json"), append=F)
  stop(displayed)
}

### Run Normalization algorithms
if (std_method_name == "log"){ # default []
  data.out = sign(data.parsed) * log2(1 + abs(data.parsed)) # Take into account normalized data with negative values
} else  if (std_method_name == "voom"){ # Default []
  require(limma)
  data.out = voom(counts=data.parsed, normalize.method="quantile", plot=F)$E
} else if (std_method_name == "cpm"){ # default []
  data.out <- 10^6 * t(data.parsed) / data.sum
  to.transpose = F
} else if (std_method_name == "tpm"){ # default []
  data.out <- t(apply(data.parsed, 2, function(x) x / data.gene.length) * 10^6) / data.sum
  to.transpose = F
} else if (std_method_name == "rpkm"){ # default [""]
  data.out <- 10^9 * t(apply(data.parsed, 2, function(x) x / data.gene.length)) / data.sum
  to.transpose = F
} else if (std_method_name == "deseq"){ # default []
  require(DESeq2)
  data.colData = data.frame(row.names=colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design=~1)
  tryCatch({
    data.dds <<- estimateSizeFactors(data.dds)
  }, error = function(err) {
    if(grepl("every gene contains at least one zero", err$message)) error.json("Every gene contains at least one zero, cannot compute log geometric means for estimating size factors.")
    error.json(err$message)
  })
  data.out <- counts(data.dds, normalized=TRUE)
  data.out <- sign(data.out) * log2(1 + abs(data.out))
} else if (std_method_name == "ruvg"){
  # Packages
  require(RUVSeq)
  
  # ERCC indexes
  data.parsed = rbind(data.parsed, data.ercc)
  data.ercc = (nrow(data.parsed) - nrow(data.ercc) + 1):nrow(data.parsed)

  # Parameters
  is_count_table = args[5]
  if(!is.null(is_count_table) & !is.na(is_count_table)) is_count_table = (is_count_table == "true")
  else is_count_table = T
  
  k = args[6]
  if(!is.null(k) & !is.na(k)) k = as.numeric(k)
  else k = 1
  
  # Run RUVg
  data.out <- RUVg(x = as.matrix(data.parsed), cIdx = data.ercc, k = k, isLog = !is_count_table)$normalizedCounts

  # Log2
  data.out <- sign(data.out) * log2(1 + abs(data.out))

  # Remove ERCC
  data.out <- data.out[-data.ercc,]
} else if(std_method_name == "seurat"){
  # Packages
  require(Seurat) # v3
  require(Matrix) # Because needed for Seurat Loom, but not embedded in Seurat v3
  
  # Parameters
  do_scale_param = args[5]
  if(!is.null(do_scale_param) & !is.na(do_scale_param)) do_scale_param = (do_scale_param == "true")
  else do_scale_param = T
  
  do_center_param = args[6]
  if(!is.null(do_center_param) & !is.na(do_center_param)) do_center_param = (do_center_param == "true")
  else do_center_param = T
  
  scale.factor_param = args[7]
  if(!is.null(scale.factor_param) & !is.na(scale.factor_param)) scale.factor_param = as.numeric(scale.factor_param)
  else scale.factor_param = 1e4
  
  scale.max_param = args[8]
  if(!is.null(scale.max_param) & !is.na(scale.max_param)) scale.max_param = as.numeric(scale.max_param)
  else scale.max_param = 10
  
  # Open the copy in writing
  data.loom <- open_with_lock(paste0(output_dir, "/output.loom"), "r+")
  
  if(do_scale_param | do_center_param){
    layer_norm = "layers/tmp"
    if(data.loom$exists(layer_norm)) data.loom$link_delete(layer_norm)
    layer_scale = output_matrix_dataset
  } else {
    layer_norm = output_matrix_dataset
  }
  
  NormalizeData.loom(object = data.loom, name = layer_norm, chunk.dims = chunk.dims_param, scale.factor = scale.factor_param, overwrite = T, display.progress = T, dataset.use = "matrix")
  if(do_scale_param | do_center_param) {
    ScaleData.loom(object = data.loom, name = layer_scale, chunk.dims = chunk.dims_param, do.scale = do_scale_param, do.center = do_center_param, scale.max = scale.max_param, normalized.data = layer_norm, overwrite = T, display.progress = T)
    data.loom$link_delete(layer_norm)
  }
  
  # Get Stats info
  nber_genes <- data.loom[[output_matrix_dataset]]$dims[2]
  nber_cells <- data.loom[[output_matrix_dataset]]$dims[1]
  nber_zeros <- sum(data.loom[[output_matrix_dataset]][,] == 0, na.rm = T)
  
  # Close the Loom file
  data.loom$close_all()
} else error.json("This normalization method is not implemented")

# Seurat directly modify the loom file, if not Seurat, add the layer to the data (original loom file)
if(std_method_name != "seurat"){
  data.loom <- open_with_lock(input_matrix_filename, "r+")
  if(data.loom$exists(output_matrix_dataset)) data.loom$link_delete(output_matrix_dataset) # Remove existing dimension reduction with same name
  if(to.transpose){
    data.loom[[output_matrix_dataset]] = t(data.out)
  } else {
    data.loom[[output_matrix_dataset]] = data.out
  }
  nber_genes <- data.loom[[output_matrix_dataset]]$dims[2]
  nber_cells <- data.loom[[output_matrix_dataset]]$dims[1]
  nber_zeros <- sum(data.loom[[output_matrix_dataset]][,] == 0, na.rm = T)
  data.loom$close_all()
} else { # If Seurat, I need to put back the created annotation in the original file
  # Run Java for copying metadata from 2 loom files
  system(paste0("java -jar ASAP.jar -T CopyMetaData -loomFrom ", output_dir, "/output.loom -loomTo ", input_matrix_filename, " -meta ", output_matrix_dataset))
  
  # Delete the temporary Loom
  file.remove(paste0(output_dir, "/output.loom"))
}

# Extra filtering step for TPM and RPKM of non-annotated genes
if(!is.null(toKeep) & sum(toKeep, na.rm = T) != nrow(data.parsed)){
  # Preparing JSON for filtering using the Java software
  stats = list()
  stats$kept_genes = which(toKeep) - 1
  write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir,"/to_keep.json"), append=F)

  # Run Java for filtering the rows in the main matrix and all metadata
  system(paste0("java -jar ASAP.jar -T FilterRows -loom ", paste0(output_dir, "/output.loom"), " -o ", paste0(output_dir, "/tmp"), " -m keep -row_names_file ", paste0(output_dir,"/to_keep.json")))
  
  # Replace the old Loom by the one we just created
  file.copy(from = paste0(output_dir, "/tmp/output.loom"), to = paste0(output_dir, "/output.loom"), overwrite = T, recursive = F, copy.mode = T, copy.date = F)
  # The output.json file is generated by the JAVA tool, so I just copy it back
  file.copy(from = paste0(output_dir, "/tmp/output.json"), to = paste0(output_dir, "/output.json"), overwrite = T, recursive = F, copy.mode = T, copy.date = F)
  
  # Delete the temporary directory
  unlink(paste0(output_dir, "/tmp/"), recursive = T, force = T)
} else {
  # Generate default JSON file
  stats <- list()
  stats$nber_rows = nber_genes
  stats$nber_cols = nber_cells
  stats$nber_zeros = nber_zeros
  write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
}
