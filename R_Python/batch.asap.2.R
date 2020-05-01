# 6) MNN doesn't handle covariates, I issue a warning if they are submitted. 
# 7) I put rather a lot of warnings (I used warning command), please check that
# they are displayed properly

# Options ---------------------------------------------------------------------
options(echo = T)
args <- commandArgs(trailingOnly = T)

# Libraries -------------------------------------------------------------------
suppressPackageStartupMessages(source("hdf5_lib.R"))
suppressPackageStartupMessages(require(jsonlite))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(scran)) # MNN
suppressPackageStartupMessages(require(sva)) # Combat

# Functions -------------------------------------------------------------------
#' collapseBatchVars
#' Collapses several batch variables into one. For example, limma can handle 
#' only 2 batch variables, then, if user submitted 3, the 3rd variable would be
#' paste together with the second one. First one is going to be untouched.
#' @param batchDF data frame with batch variables, rows are samples and columns
#'                are batch variables
#' @param colStart number of batch variables a method can handle.
#' @example
#' for limma (2 variables): collapseBatchVars(batchTab, 2)
#' for combat (1 variable): collapseBatchVars(batchTab, 1)
collapseBatchVars <- function(batchDF, colStart = 2) {
  result <- data.frame(batchDF)
  if (ncol(batchDF) > colStart) {
    msg <- paste("The method can't handle more than", colStart,
                 "categorical variables. Remaining categorical variables are assumed to be ", 
                 "additive & collapsed")
    warning(msg)
	data.warnings[[3]] <- list(name = msg, description = msg)
    result[, colStart] <- apply(result[, colStart:ncol(batchDF)], 1, paste, 
                                collapse = ', ')
  }
  return(result)
}

#' combatBatchCorrect
#' Corrects batch effect with use of sva (aka Combat) package
#' @param dataTab NORMALIZED expression data frame, columns are cells
#' @param batchTab data frame with batch data, default: NULL
#' @param covarsTab data frame with covariates, default: NULL
#' @return data frame with normalized data
#' @note can handle ONLY 1 BATCH EFFECT variable. UNLIMITED COVARIATES
#' @note from the function description about input table: The input data are 
#'       assumed to be cleaned and normalized before batch effect removal.
combatBatchCorrect <- function(dataTab, batchTab = NULL, covarsTab = NULL) {
  # correction with only covariates isn't possible
  if (is.null(batchTab) & !is.null(covarsTab)) {
    msg <- paste("Combat requires at least one categorical variable. Removal of",
                 "only numerical covariates are not supported by",
                 "Combat. Try using limma instead.")
    error.json(msg)
    result <- NULL
  }
  
  if (!is.null(batchTab)) {
    # Convert to vector / collapse as Combat can only handle one batch effect
    # variable. Additive effect is assumed.
    if (ncol(batchTab) > 1) {
      batchVect <- collapseBatchVars(batchTab, 1)[, 1]
    } else {
      batchVect <- batchTab[, 1]
    }
    
    # create formula for covariates, if they are submitted
    if (!is.null(covarsTab)) {
      modelAsStr <- paste("~ ", paste0(colnames(covarsTab), collapse = ' + '))
      modelEqv <- as.formula(modelAsStr)
      modelMatr <- model.matrix(modelEqv, covarsTab)
      
      result <- ComBat(as.matrix(dataTab), batchVect, modelMatr)
    } else {
      result <- ComBat(as.matrix(dataTab), batchVect)
    }
  }
  rownames(result) <- rownames(dataTab)
}

#' limmaBatchCorrect
#' Corrects batch effect with use of limma package
#' @param dataTab NORMALIZED expression data frame, columns are cells
#' @param batchTab data frame with batch data, default: NULL
#' @param covarsTab data frame with covariates, default: NULL
#' @return data frame with normalized data
#' @note can handle 2 BATCH EFFECT variables. UNLIMITED COVARIATES
#' @note from the function description about input table: numeric matrix, or 
#' any data object that can be processed by getEAWP containing log-expression 
#' values for a series of samples. Rows correspond to probes and columns to 
#' samples.
limmaBatchCorrect <- function(dataTab, batchTab = NULL, covarsTab = NULL) {
  if (is.null(batchTab) & !is.null(covarsTab)) {
    result <- removeBatchEffect(dataTab, covariates = covarsTab)
  }
  if (!is.null(batchTab)) {
    # Limma can't handle > 2 batch variables. So, if more than 2 batch variables
    # I assume additive effect and pull them together.
    if (ncol(batchTab) > 2) {
      batchTab <- collapseBatchVars(batchTab, 2)
      msg <- paste("The first batch variable (", colnames(batchTab)[1], ")",
                   "was assumed to be the most important one and isn't",
                   "collapsed")
      warning(msg)
	  data.warnings[[4]] <- list(name = msg, description = msg)
    }
    if (ncol(batchTab) > 1) {
      result <- removeBatchEffect(dataTab, batchTab[, 1], batchTab[, 2], 
                                  covarsTab)
    } else {
      result <- removeBatchEffect(dataTab, batchTab[, 1], 
                                  covariates = covarsTab)
    }
  }
  rownames(result) <- rownames(dataTab)
  result
}

#' mnnBatchCorrect
#' Corrects batch effect with use of scran (aka MNN) package
#' @param dataTab NORMALIZED expression data frame, columns are cells
#' @param batchTab data frame with batch data, default: NULL
#' @param covarsTab data frame with covariates, default: NULL
#' @param k integer scalar specifying the number of nearest neighbors to 
#'          consider when identifying mutual nearest neighbors.
#' @param sigma numeric scalar specifying the bandwidth of the Gaussian 
#'              smoothing kernel used to compute the correction vector for each
#'              cell.
#' @return data frame with normalized data
#' @note can handle 1 BATCH EFFECT variables. NO COVARIATES.
#' @note from the function description about input table: The input expression
#' values should generally be log-transformed, e.g., log-counts, see normalize
#' for details. They should also be normalized within each data set to remove 
#' cell-specific biases in capture efficiency and sequencing depth. 
mnnBatchCorrect <- function(dataTab, batchTab = NULL, covarsTab = NULL, 
                            mnnK = 20, mnnSigma = 0.1) {
  # correction with only covariates isn't possible
  if (!is.null(covarsTab)) error.json("Removal of numerical covariates are not supported by MNN. Try using limma instead.")
  if (is.null(batchTab)) error.json("MNN works only if you supply at least one categorical covariate!")
  
  # Convert to vector / collapse as MNN can only handle one batch effect 
  # variable. Additive effect is assumed.
  if (ncol(batchTab) > 1) {
	batchVect <- as.factor(collapseBatchVars(batchTab, 1)[, 1])
  } else {
    batchVect <- as.factor(batchTab[, 1])
  }

  # MNN wants banch of matrices as input
  batchMatrx <- lapply(levels(batchVect), function(x) dataTab[, which(batchVect == x)])
  batchMatrx <- lapply(batchMatrx, as.matrix)
  result <- do.call(mnnCorrect, c(batchMatrx, list(k = mnnK, sigma = mnnSigma)))
  result <- result$corrected
	
  # recover column names
  result <- do.call(cbind, result)
  colnames(result) <- do.call(c, lapply(batchMatrx, colnames))

  # recover order
  result <- result[, colnames(dataTab)]

  rownames(result) <- rownames(dataTab)
  result
}

# Default Parameters  ---------------------------------------------------------
set.seed(42)
# names of implemented methods
batchRemoveMethods <- c('limma', 'combat', 'mnn') 

input_matrix_filename <- args[1]
output_dir <- args[2]
std_method_name <- args[3]
input_matrix_dataset <- args[4]
output_matrix_dataset <- args[5]
batch_dataset <- args[6]
data.warnings <- NULL
time_idle <- 0

if (std_method_name == 'mnn') {
  mnn_k <- 20
  if (!is.na(args[7])) {
    mnn_k <- as.integer(args[7])
  } else {
	data.warnings[[1]] <- list(name = "No k for MNN was submitted, using default one (20)", description = "No k for MNN was submitted, using default one (20)")
  }
  if (!is.na(args[8])) {
    mnn_sigma <- as.numeric(args[8])
  } else {
    data.warnings[[2]] <- list(name = "No sigma for MNN was submitted, using default one (0.1)", description = "No sigma for MNN was submitted, using default one (0.1)")
  }
}

# check, that method is implemented
if (!std_method_name %in% batchRemoveMethods) {
  error.json(paste0(std_method_name, " method is not implemented"))
}

# check that batch variables are submitted
if(is.null(batch_dataset) || is.na(batch_dataset) || batch_dataset == "" ||
   batch_dataset == "null"){ 
  error.json("To remove batch effect you need to submit batch variables!")
} else {
  batch_dataset <- unique(strsplit(batch_dataset, split = ",")[[1]])
  message(length(batch_dataset), " covariates detected. Will add covariates to the model.")
}

# Read-in LOOM file -----------------------------------------------------------
# Open the existing Loom in read-only mode and recuperate the infos 
# !!! not optimized for OUT-OF-RAM computation !!!
data.loom <- open_with_lock(input_matrix_filename, "r")
# parsed data
data.parsed <- fetch_dataset(data.loom, input_matrix_dataset, transpose = T)
if(is.null(data.parsed)) {
  error.json("This dataset does not exist in the Loom file")
}
# gene names
data.ens <- fetch_dataset(data.loom, "/row_attrs/Accession")
data.gene <- fetch_dataset(data.loom, "/row_attrs/Gene")

# Parse batch effect vars: remove duplicates and pass numeric as covariates----
# remove non-unique batch variables
data.batch <- lapply(batch_dataset, function(x) fetch_dataset(data.loom, x))
close_file(data.loom)

# separate pass numeric batch effect as covariates
batchClass <- sapply(data.batch, class)
numericBatch <- !batchClass %in% c('factor', 'logical', 'character')
if (any(numericBatch)) {
  msg <- paste("Detected numeric batch effect variable(s) ",
               paste0(batch_dataset[numericBatch], collapse = ','), 
               ". Will treat them as covariates!")
  warning(msg)
  data.covars <- data.batch[!batchClass %in% 
                              c('factor', 'logical', 'character')]
  data.batch <- data.batch[batchClass %in% c('factor', 'logical', 'character')]
}

# convert to data.frame
if (length(data.batch) != 0) {
  data.batch <- as.data.frame(do.call(cbind, data.batch))
  colnames(data.batch) <- batch_dataset[!numericBatch]
  # replace non-alphanumeric characters from the column names, so we could use
  # names in the formula
  colnames(data.batch) <- gsub('.*/', '', colnames(data.batch))
  colnames(data.batch) <- gsub('^_', '', colnames(data.batch))
  message(paste("Batch effect from: ",
                paste(colnames(data.batch), collapse = ', '),
                "will be removed"))
} else {
  data.batch <- NULL
}
if (exists("data.covars")) {
  data.covars <- as.data.frame(do.call(cbind, data.covars))
  colnames(data.covars) <- batch_dataset[numericBatch]
  # replace non-alphanumeric characters from the column names, so we could use
  # names in the formula
  colnames(data.covars) <- gsub('.*/', '', colnames(data.covars))
  colnames(data.covars) <- gsub('^_', '', colnames(data.covars))
  message(paste("Effect from: ", paste(colnames(data.covars), collapse = ', '),
                "will be regressed out as covariate"))
} else {
  data.covars <- NULL
}

# fool defense
if (is.null(data.batch) & is.null(data.covars)) {
  msg <- paste("No batch or covariates to regress out is detected.",
               "Something went very wrong")
  error.json(msg)
}

# Perform batch correction ----------------------------------------------------
data.out <- switch(std_method_name, 
                   "limma" = limmaBatchCorrect(data.parsed, data.batch,
                                               data.covars),
                   "combat" = combatBatchCorrect(data.parsed, data.batch, 
                                                 data.covars),
                   "mnn" = mnnBatchCorrect(data.parsed, data.batch, 
                                           data.covars, mnn_k, mnn_sigma))

# Output results --------------------------------------------------------------
# Open Loom in writing mode for writing results
data.loom <- open_with_lock(input_matrix_filename, "r+")
add_matrix_dataset(handle = data.loom, dataset_path = output_matrix_dataset, dataset_object = t(data.out))
close_all()

# Generate default JSON file
stats <- list()
stats$nber_rows = nrow(data.out)
stats$nber_cols = ncol(data.out)
if(!is.null(data.warnings)) stats$warnings = list(data.warnings)
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)