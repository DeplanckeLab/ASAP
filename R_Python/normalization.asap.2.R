### ASAP Normalization script
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
output_matrix_dataset <- args[4]
toKeep <- NULL
time_idle <- 0

#### Functions ####
LogNorm <- function(data, scale_factor, display_progress = TRUE) {
  .Call('_Seurat_LogNorm', PACKAGE = 'Seurat', data, scale_factor, display_progress)
}

FastRowScale <- function(mat, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastRowScale', PACKAGE = 'Seurat', mat, scale, center, scale_max, display_progress)
}

# Read Loom
data.loom <- open_with_lock(input_matrix_filename, "r") # Block access to the file (nobody can write if we copy)
data.parsed <- fetch_dataset(data.loom, "/matrix", transpose = T) # t() because loomR returns the t() of the correct matrix we want
if(std_method_name %in% c("tpm", "rpkm")) {
  data.gene.length <- fetch_dataset(data.loom, "/row_attrs/_SumExonLength")
  toKeep <- data.gene.length > 0
  data.gene.length[!toKeep] = NA
}
if(std_method_name %in% c("cpm", "tpm", "rpkm")) data.sum <- fetch_dataset(data.loom, "/col_attrs/_Depth")
if(std_method_name == "ruvg") {
  data.ercc <- fetch_dataset(data.loom, "/col_attrs/_ERCCs")
  if(is.null(data.ercc) || nrow(data.ercc) == 0) error.json(displayed = "RUVg method can only be ran with ERCCs (spike-ins) used as 'negative control genes'")
}
close_all()

### Run Normalization algorithms
if (std_method_name == "log"){ # default []
  data.out = t(sign(data.parsed) * log2(1 + abs(data.parsed))) # Take into account normalized data with negative values
} else  if (std_method_name == "voom"){ # Default []
  require(limma)
  data.out = t(voom(counts=data.parsed, normalize.method="quantile", plot=F)$E)
} else if (std_method_name == "cpm"){ # default []
  data.out <- 10^6 * t(data.parsed) / data.sum
  data.out[is.na(data.out)] <- 0 # If data.sum = 0 it means that the full col is 0 no?
} else if (std_method_name == "tpm"){ # default []
  data.out <- t(apply(data.parsed, 2, function(x) x / data.gene.length) * 10^6) / data.sum
  data.out[is.na(data.out)] <- 0 # If data.sum = 0 it means that the full col is 0 no?
} else if (std_method_name == "rpkm"){ # default [""]
  data.out <- 10^9 * t(apply(data.parsed, 2, function(x) x / data.gene.length)) / data.sum
  data.out[is.na(data.out)] <- 0 # If data.sum = 0 it means that the full col is 0 no?
} else if (std_method_name == "deseq"){ # default []
  require(DESeq2)
  data.colData = data.frame(row.names=colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design=~1)
  tryCatch({
    data.dds <<- estimateSizeFactors(data.dds)
  }, error = function(err) {
    if(grepl("every gene contains at least one zero", err$message)) error.json("Every gene contains at least one zero, cannot compute log geometric means for estimating size factors. Maybe filter more your dataset to remove low expressed genes/weak cells.")
    error.json(err$message)
  })
  data.out <- counts(data.dds, normalized=TRUE)
  data.out <- t(sign(data.out) * log2(1 + abs(data.out)))
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
  data.out <- t(data.out[-data.ercc,])
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

  data.out <- as.matrix(LogNorm(data = Matrix(data = t(x = data.parsed), sparse = T), scale_factor = scale.factor_param, display_progress = F))
  if(do_scale_param | do_center_param) data.out <- FastRowScale(mat = data.out, scale = do_scale_param, center = do_center_param, scale_max = scale.max_param, display_progress = F)
} else error.json("This normalization method is not implemented")

# Writing the dataset in the loom file
data.loom <- open_with_lock(input_matrix_filename, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_matrix_dataset, dataset_object = data.out)
close_all()

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
  stats$nber_rows = ncol(data.out)
  stats$nber_cols = nrow(data.out)
  write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
}
