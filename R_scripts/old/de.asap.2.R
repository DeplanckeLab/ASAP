### ASAP DE script

# Libraries -------------------------------------------------------------------
source("hdf5_lib.R")
require(jsonlite)
require(ggplot2)
require(plotly)
require(stats) # [MASHA]

# Functions -------------------------------------------------------------------
serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty = T, digits = 3)
}

# [MASHA]
#' bimodLikData
#' Copied from Seurat
#' internal function to run mcdavid et al. DE test
#' @importFrom stats sd dnorm
bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- MinMax(data = length(x = x2) / length(x = x), min = 1e-5, max = (1 - 1e-5))
  likA <- length(x1) * log(1 - xal)
  if (length(x2) < 2) {
    mysd <- 1
  } else {
    mysd <- sd(x2)
  }
  likB <- length(x2) * log(xal) +
          sum(dnorm(x2, mean = mean(x2), sd = mysd, log = T))
  return(likA + likB)
}

# [MASHA]
#' DifferentialLRT
#' Copied from Seurat
#' internal function to run mcdavid et al. DE test
#' @importFrom stats pchisq
DifferentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimodLikData(x)
  lrtY <- bimodLikData(y)
  lrtZ <- bimodLikData(c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

# [MASHA]
#' applyNegbinomOrPoisson 
#' Applies negative binomial or poisson test to expression. Adapted from Seurat
#' @param testName name of the test
#' @param geneName character, gene name
#' @param minPct double, minimal number of cells in the group
#' @param exprVals double, expression
#' @param testMetaData data.frame with 1st column called group, 
#'                     and rest - covariates
#' @return p-value
applyNegbinomOrPoisson <- function(testName, geneName, minPct, exprVals, 
                                   testMetaData) {
  testData <- data.frame(Expr = unlist(exprVals), testMetaData)

  # check that gene is expressed in specified number of cells in one group
  gr1n <- sum(testData$Expr[testData$group == levels(testData$group)[1]] > 0)
  gr1n <- gr1n / nrow(testData)
  gr2n <- sum(testData$Expr[testData$group == levels(testData$group)[2]] > 0)
  gr2n <- gr2n / nrow(testData)
  if (gr1n < minPct && gr2n < minPct) {
    warning(paste0( "Skipping gene --- ", geneName,  ". Fewer than ", minPct,
                    " cells in both clusters."))
    return(NA)
  }
  # check that variance between groups is not 0
  if (var(testData$Expr) == 0) {
    warning(paste0("Skipping gene -- ", geneName, 
                   ". No variance in expression between the two clusters."))
    return(NA)
  }
  # create formula:
  modelAsStr <- "Expr ~ group"
  if (ncol(testData) > 2) {
    modelAsStr <- paste(modelAsStr, " + ", 
                        paste(colnames(testData)[-2:-1], collapse = ' + '))
  }
  modelEqv <- formula(modelAsStr)
  
  # perform test
  if (testName == "negbinom-seurat") {
    try(p.est <- summary(glm.nb(formula = modelEqv, 
                                data = testData))$coef[2, 4],
        silent = TRUE)
    return(p.est)
  }
  if (testName == "poisson-seurat") {
    try(p.est <- summary(glm(formula = modelEqv, family = "poisson",
                             data = testData))$coef[2, 4],
        silent = TRUE)
    return(p.est)
  }
  if (!testName %in% c("negbinom-seurat", "poisson-seurat")) {
    return(NA)
  }
}

#  [MASHA]
#' Imported from Seurat
#' Calculate the mean of logged values
#' Calculate mean of logged values in non-log space (return answer in 
#' log-space)
#' @param x A vector of values
#' @param ... Other arguments (not used)
#' @return Returns the mean in log-space
#' @export
#' @examples
#' ExpMean(x = c(1, 2, 3))
ExpMean <- function(x, ...) {
  x <- as.numeric(x)
  if (inherits(x = x, what = 'AnyMatrix')) {
    return(apply(x, 1, function(i) log(x = mean(x = exp(x = i) - 1) + 1)))
  } else {
    return(log(mean(exp(x) - 1) + 1))
  }
}

#  [MASHA]
#' Imported from Seurat
#' internal function to calculate AUC values
#' @importFrom ROCR prediction performance
#'
DifferentialAUC <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  prediction.use <- prediction(predictions = c(x, y),
                               labels = c(rep(x = 1, length(x = x)), 
                                          rep(x = 0, length(x = y))),
                               label.ordering = 0:1)
  perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}

#  [MASHA]
#' Imported from Seurat
#' internal function to calculate AUC values
#' @importFrom pbapply pblapply
AUCMarkerTest <- function(data1, data2, mygenes) {
  myAUC <- sapply(mygenes, function(x) DifferentialAUC(data1[x, ], data2[x, ]))
  myAUC[is.na(myAUC)] <- 0
  avg_diff <- sapply(mygenes, 
                     function(x) ExpMean(data1[x, ]) - ExpMean(data2[x, ]))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  toRet
}

# Options ---------------------------------------------------------------------
# Second option is for future_sapply (100Gb allowed here)
options(echo = T, future.globals.maxSize = 100E9) 
args <- commandArgs(trailingOnly = T)

# Default Parameters  ---------------------------------------------------------
set.seed(42)
time_idle <- 0
# [MASHA] list of implemented DE methods
availableMethods <- c('limma', 'deseq2', 'wilcox_seurat', 'ttest_seurat', 
                      'bimod_seurat', 'negbinom_seurat', 'poisson_seurat', 
                      'mast_seurat', 'lr_seurat') 

input_matrix_filename <- args[1]
output_dir <- args[2]
std_method_name <- args[3]
input_matrix_dataset <- args[4]
output_matrix_dataset <- args[5]
batch_dataset <- args[6]
group_dataset <- args[7]
group_1 <- args[8]
group_2 <- args[9]
is_count_table <- args[10]
min_pct <- as.numeric(args[11])
min_diff_pct <- args[12]
fc_threshold <- as.numeric(args[13])
max_cells_per_ident <- args[14]
nb.cores <- args[15] # DESEQ2 specific # Note: if length(args) < 15, args[15] = NA

# [MASHA] test params
#input_matrix_filename <- "/data/vincent/epmt3d_de_35963_output.loom" #args[1]
#output_dir <- "/data/vincent" #args[2]
#std_method_name <- "wilcox_seurat" #args[3]
#input_matrix_dataset <- "/layers/norm_1_asap_seurat" #args[4]
#output_matrix_dataset <- "/row_attrs/_de_97_wilcox_seurat" #args[5]
# [MASHA] I had to vectorize here, otherwise I couldn't extract data from loom
## via fetch_dataset
#batch_dataset <- "/col_attrs/_StableID" #args[6]
#group_dataset <- "/col_attrs/_clust_2_seurat" # args[7]
#group_1 <- "4" #args[8]
#group_2 <- "null" #args[9]
#is_count_table <- "false" #args[10]
#min_pct <- 0.1  # as.numeric(args[11])
#min_diff_pct <- "null" # as.numeric(args[12])
#fc_threshold <- 1.3 # as.numeric(args[13])
#max_cells_per_ident <- "null" #args[14]
#nb.cores <- NA # args[15] # DESEQ2 specific

# [MASHA]
if (!std_method_name %in% availableMethods) {
  error.json(paste0(std_method_name, " method is not implemented"))
}

if(group_1 == group_2) {
  error.json("Cannot compute DE from the same group.")
}

if(is_count_table == "true") { 
  is_count_table <- T
} else if(is_count_table == "false") {
  is_count_table <- F
} else {
  error.json("is_count_table should be 'true' or 'false'")
}

if(is.null(min_diff_pct) || is.na(min_diff_pct) || min_diff_pct == "" || min_diff_pct == "null"){ 
  min_diff_pct <- -Inf
} else {
  min_diff_pct <- as.numeric(min_diff_pct)
}

if(is.null(max_cells_per_ident) || is.na(max_cells_per_ident) || max_cells_per_ident == "" || max_cells_per_ident == "null") { 
  max_cells_per_ident <- Inf
} else {
  max_cells_per_ident <- as.numeric(max_cells_per_ident)
}
data.warnings <- NULL

# Read-in LOOM file -----------------------------------------------------------
# Open the existing Loom in read-only mode and recuperate the infos 
# !!! not optimized for OUT-OF-RAM computation !!!
data.loom <- open_with_lock(input_matrix_filename, "r")
data.parsed <- fetch_dataset(data.loom, input_matrix_dataset, transpose = T)
data.groups <- fetch_dataset(data.loom, group_dataset)
data.ens <- fetch_dataset(data.loom, "/row_attrs/Accession")
data.gene <- fetch_dataset(data.loom, "/row_attrs/Gene")
if(is.null(batch_dataset) || is.na(batch_dataset) || batch_dataset == "" || batch_dataset == "null"){ 
  data.batch <- NULL
} else {
  batch_dataset <- unique(strsplit(batch_dataset, split = ",")[[1]])
  message(length(batch_dataset), " covariates detected. Will add covariates to the model.")
  data.batch <- lapply(batch_dataset, function(x) fetch_dataset(data.loom, x))
  data.batch <- lapply(data.batch, as.data.frame)
  data.batch <- do.call(cbind, data.batch)
  # replace non-alphanumeric characters from the column names, so we could use
  # names in the formula
  colnames(data.batch) <- gsub('.*/', '', batch_dataset)
  colnames(data.batch) <- gsub('^_', '', colnames(data.batch))
  message(paste("Covariates added to the model:", paste(colnames(data.batch), collapse = ', ')))
}
close_file(data.loom)

if(is.null(data.parsed)) {
  error.json("This dataset does not exist in the Loom file")
}
if(is.null(data.groups)) {
  error.json(paste0(group_dataset, " does not exist in the Loom file"))
}

# Are data log2 or not? -------------------------------------------------------
# Assessing log2ity of the data (could also be ln or log10 that should still be ok)
is_log2 <- T
compute_FC <- T

# Arbitrary cutoff, could probably be set lower than that
if (is_count_table | max(data.parsed) > 30) {
  is_log2 <- F
}
if(min(data.parsed) < 0) {
  compute_FC <- F
  is_log2 <- F
  is_count_table <- F
  msg <- "Data has negative values, no FC will be computed"
  descrMsg <- "Select another dataset if you want to compute FC."
  data.warnings[[1]] <- list(name = msg, description = descrMsg)
}

# Handle groups ---------------------------------------------------------------
groups <- unique(data.groups)
data.group <- rep(NA, ncol(data.parsed))
data.group[data.groups == group_1] <- 1
if(group_2 != "null"){
  msg <- paste("Two group files found.", "Performing DE [G1 vs G2] and discarding other samples")
  message(msg)
  data.group[data.groups == group_2] <- 2
} else {
  msg <- paste("One group file found.", "Performing this group against all other samples [Marker Genes]")
  message(msg)
  data.group[data.groups != group_1] <- 2
}

# Filter NA groups // keep only cells in the designated groups
toKeep <- which(!is.na(data.group))
data.parsed <- data.parsed[, toKeep]
data.batch <- data.batch[toKeep, , drop=F]
data.group <- data.group[toKeep]
cells_1 <- data.group == 1
cells_2 <- data.group == 2
features <- 1:nrow(data.parsed)

# Check if computable
if(sum(data.group == 1) < 3) {
  error.json("Group 1 should contain at least 3 samples")
}
if(sum(data.group == 2) < 3) {
  error.json("Group 2 should contain at least 3 samples")
}

if(!is.null(nb.cores) & !is.na(nb.cores) & nb.cores != "" & nb.cores != "null") {
  nb.cores <- as.numeric(nb.cores)
} else {
  nb.cores <- 8
}

#[Seurat] Feature selection (based on percentages) to speed up computation ----
if(is_count_table | is_log2){ # Because checking genes > 0 as "detected"
  pct.1 <- rowSums(x = data.parsed[, cells_1, drop = FALSE] > 0)
  pct.1 <- round(x = pct.1 / sum(x = cells_1), digits = 3)
  pct.2 <- rowSums(x = data.parsed[, cells_2, drop = FALSE] > 0)
  pct.2 <- round(x = pct.2 / sum(x = cells_2), digits = 3)
  
  data.alpha <- cbind(pct.1, pct.2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  features <- which(x = alpha.min > min_pct)
  if(length(x = features) == 0) {
    error.json("No features pass min_pct threshold")
  }
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
  features <- which(x = alpha.min > min_pct & alpha.diff > min_diff_pct)
  if(length(x = features) == 0) {
    error.json("No features pass min_diff_pct threshold")
  }
} else {
  data.warnings[[2]] <- list(name = "Data is neither a count matrix nor detected as logged. No pre-filtering (using &#8217;Min Pct&#8217; and &#8217;Min Diff Pct&#8217;) was performed", 
                             description = "Select another dataset if you want to perform pre-filtering (because looking for &#8217;detected genes&#8217; i.e. with values of expression &gt; 0).")
}

#[Seurat] Gene selection (on average difference/FC) to speed up computation----
total.diff <- NULL
if(compute_FC){
  # If also !is_count_table then I don't know the data type, so I still log
  if(!is_log2){ 
    data.1 <- apply(X = data.parsed[, cells_1, drop = FALSE], MARGIN = 1, 
                    FUN = function(x) { m = mean(x); 
                                        return(sign(m) * log2(abs(m) + 1))})
    data.2 <- apply(X = data.parsed[, cells_2, drop = FALSE], MARGIN = 1,
                    FUN = function(x) { m = mean(x); 
                                        return(sign(m) * log2(abs(m) + 1))})
  } else { # should contain negative values
    data.1 <- apply(X = data.parsed[, cells_1, drop = FALSE], MARGIN = 1, 
                    FUN = function(x) return(log2(mean(2^x - 1) + 1)))
    data.2 <- apply(X = data.parsed[, cells_2, drop = FALSE], MARGIN = 1, 
                    FUN = function(x) return(log2(mean(2^x - 1) + 1)))
  } 
  total.diff <- data.1 - data.2
  features.diff <- which(abs(total.diff) > log2(fc_threshold))
  features <- intersect(x = features, y = features.diff)
  if(length(x = features) == 0) {
    error.json("No features pass 'FC threshold'")
  }
}

#[Seurat] Downsampling cells if too many per group ----------------------------
if (max_cells_per_ident < Inf) {
  if (sum(cells_1) > max_cells_per_ident) {
    cells_1_kept <- sample(x = which(cells_1), size = max_cells_per_ident)
  } else { 
    cells_1_kept = which(cells_1)
  }
  if (sum(cells_2) > max_cells_per_ident) {
    cells_2_kept <- sample(x = which(cells_2), size = max_cells_per_ident)
  } else {
    cells_2_kept = which(cells_2)
  }
  
  # Filter the datasets
  toKeep = sort(c(cells_1_kept, cells_2_kept))
  data.parsed <- data.parsed[, toKeep] # We keep all features
  data.batch <- data.batch[toKeep, ,drop=F] # [MASHA]
  data.group <- data.group[toKeep]
  cells_1 <- data.group == 1
  cells_2 <- data.group == 2
}

# Run LIMMA -------------------------------------------------------------------
if(std_method_name == "limma"){
  # Package
  require(limma)
  
  # Create design matrix from groups
  if(is.null(data.batch)){
    data.design <- model.matrix(~ 0 + data.group)
  } else { # [MASHA] incorporate the covariats into the model
    # model matrix of interest
    data.colData <- data.frame(group = factor(data.group), data.batch)
    # create formula
    modelAsStr <- paste("~ 0 + ", paste0(colnames(data.colData), collapse = ' + '))
    data.design <- model.matrix(formula(modelAsStr), data.colData)
    # define contrast matrix, COI = comparison of interest
    contrastMatr <- makeContrasts(COI = paste(colnames(data.design)[1:2], collapse = '-'), levels = data.design)
  }
  
  # Test if the data is count or normalized already
  if(is_count_table){ # if data is a count matrix
    message("Input data is a count table. Voom Normalizing...")
    data.voom <- voom(counts = data.parsed[features, ], normalize.method = "quantile", plot = F)$E
    fit <- NULL
    tryCatch({ fit <- lmFit(data.voom, data.design)},
             warning = function(war) {
               if(grepl("Partial NA coefficients", war$message)) {
                 msg <- paste("Your batches and selections are confounded.",
                              "There is no way to properly analyze the",
                              "experiment in that case.")
                 error.json(msg)
               }
             })
  } else { 
    fit <- NULL
    tryCatch({ if(!is_log2 | !compute_FC) {
        fit <- lmFit(data.parsed[features, ], data.design)
        } else {
          fit <- lmFit(log2(1 + data.parsed[features,]), data.design)
      }
    }, warning = function(war) {
      if(grepl("Partial NA coefficients", war$message)) {
        msg <- paste("Your batches and selections are confounded.",
                     "There is no way to properly analyze the",
                     "experiment in that case.")
        error.json(msg)
      }
    })
  }
  
  # [MASHA] Run the DE model with account for covariats
  if (is.null(data.batch)) {
    fit  <- eBayes(fit, trend = TRUE)
  } else {
    fit  <- contrasts.fit(fit, contrastMatr)
    fit  <- eBayes(fit, trend = TRUE)
  }

  # Create the table with all DE results
  DE <- topTable(fit, n = Inf, adjust = "fdr", sort.by = "none")
  
  # Create output
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) <- c("logFC", "pval", "FDR", "AveG1", "AveG2")
  # [MASHA] I put row.names(DE) for safety
  data.out[as.integer(row.names(DE)), "logFC"] <- DE$logFC 
  data.out[as.integer(row.names(DE)), "pval"] <- DE$P.Value
  data.out[as.integer(row.names(DE)), "FDR"] <- DE$adj.P.Val
  data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
}

# Run DESEQ2 ------------------------------------------------------------------ 
if (std_method_name == "deseq2"){ # Assuming data is a count table
  # Packages
  require(DESeq2)
  require(BiocParallel)
  
  # Check
  if(!is_count_table) {
    error.json("Data should be a count table in order to run DESeq2")
  }
  
  # create DDS
  data.dds <- NULL
  data.colData <- data.frame(group = factor(data.group))
  if(is.null(data.batch)){
    modelEqv <- formula('~ group')
  } else { # [MASHA] added covariats to the model
    # model matrix, variable of interest at the end as vignette says
    data.colData <- cbind(data.batch, data.colData)
    # create formula
    modelAsStr <- paste("~ ", paste0(colnames(data.colData), collapse = ' + '))
    modelEqv <- formula(modelAsStr)
    # contrasts of interest - for future
    COI <- c('group', levels(data.colData$group)[1], 
             levels(data.colData$group)[2])
  }
  tryCatch({
    data.dds <- DESeqDataSetFromMatrix(data.parsed[features, ],
                                       colData = data.colData, 
                                       design = modelEqv)
  }, error = function(err) {
    if(grepl("the model matrix is not full rank", err$message)) {
      msg <- paste0("Seems that your batches and selections are confounded.",
                    "There is no way to properly analyze the experiment in",
                    "that case.")
      error.json(msg)
    }
    error.json(err$message)
  })
  
  tryCatch({
    if(nb.cores > 1) {
      register(MulticoreParam(nb.cores))
      data.dds <- DESeq(data.dds, parallel = T)
    } else {
      data.dds <- DESeq(data.dds, parallel = F)
    }
  }, error = function(err) {
    if(grepl("every gene contains at least one zero", err$message)) {
      msg <- paste("Every gene contains at least one zero, cannot compute log",
                   "geometric means for estimating size factors.")
      error.json(msg)
    }
    if(!is.null(err)) {
      error.json(err$message)
    }
  })
  
  # Obtain the DE results
  if (is.null(data.batch)) {
    DE <- results(data.dds)
  } else { # [MASHA]
    DE <- results(data.dds, contrast = COI)
  }
  
  # Create output
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  # [MASHA] I put row.names(DE) for safety
  data.out[as.integer(row.names(DE)), "logFC"] <- DE$log2FoldChange
  data.out[as.integer(row.names(DE)), "pval"] <- DE$pvalue
  data.out[as.integer(row.names(DE)), "FDR"] <- DE$padj
  data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
}

# [Seurat] Run wilcox ---------------------------------------------------------
if(std_method_name == "wilcox_seurat"){ # Seurat function
  require(Seurat)
  require(future.apply)
  
  # parallelization
  plan(multiprocess, workers = nb.cores) # parallelize
  
  # Design matrix
  data.colData <- data.frame(GOI = factor(data.group))
  # Perform the test
  if(is.null(data.batch)) { # without covariates
    testPvals <- future_sapply(X = features, 
                               FUN = function(x) {
                                 geneExpr <- unlist(data.parsed[x, ])
                                 res <- wilcox.test(geneExpr ~ data.colData$GOI)
                                 res <- res$p.value
                                 return(res)
                               })
  } else { # [MASHA] Wilcoxon test with covariats
    require(sanon)
	  data.warnings[[3]] <- list(name = "Wilcoxon test with covariates isn't default for Seurat: package 'sanon' was used", description = "Wilcoxon test with covariates is much slower")
    
    data.colData <- data.frame(data.colData, data.batch)
    
    # create formula:
    # 1) add covariates
    modelAsStr <- paste0("Expr ~ covar(", colnames(data.batch), ")", collapse = " + ")
    
    # 2) add group of interest
    modelAsStr <- paste0(modelAsStr, ' + grp(GOI)')
    modelEqv <- formula(modelAsStr)
    testPvals <- future_sapply(X = features, 
                               FUN = function(x) {
                                 testData <- data.frame(Expr = unlist(data.parsed[x, ]), data.colData)
                                 res <- summary(sanon(modelEqv, data = testData))
                                 res <-res$coefficients[, 4]
                               })
  }
  
  # Create output table
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features, "pval"] <- testPvals
  data.out[features, "FDR"] <- p.adjust(p = data.out[features, "pval"], method = "fdr")
  if(compute_FC) data.out[, "logFC"] <- total.diff
  data.out[, "AveG1"] = rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] = rowMeans(data.parsed[, cells_2, drop = F])
} 

# [Seurat] Run t-test by [MASHA] ----------------------------------------------
if(std_method_name == "ttest_seurat"){ # Seurat function DiffTTest as base
  require(Seurat)
  require(future.apply)
  
  # Apply the test
  data.colData <- data.frame(group = factor(data.group))
  if(is.null(data.batch)) { # t-test without covariats
    testPvals <- future_sapply(X = features,
                               function(x) t.test(unlist(data.parsed[x, ]) ~ 
                                                  data.colData$group)$p.value)
  } else { 
    data.colData <- data.frame(data.colData, data.batch)
    # create formula:
    modelAsStr <- paste("Expr ~ ", paste0(colnames(data.colData), 
                                          collapse = ' + '))
    modelEqv <- formula(modelAsStr)
    # lm and t.test are equal
    testPvals <- future_sapply(X = features,
                               function(x) {
                               testData <- data.frame(Expr = unlist(data.parsed[x, ]),
                                                      data.colData)
                               summary(lm(modelEqv, testData))$coefficients[2, 4]
                               })
  }
  
  # Create output template
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) <- c("logFC", "pval", "FDR", "AveG1", "AveG2")    
  data.out[features, "pval"] <- testPvals
  data.out[features, "FDR"] <- p.adjust(p = data.out[features, "pval"],
                                        method = "fdr")
  if(compute_FC) data.out[, "logFC"] <- total.diff
  data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
} 

# [Seurat] Run bimod by [MASHA] -----------------------------------------------
if(std_method_name == "bimod_seurat") { # Seurat function DiffExpTest as a base
  require(Seurat)
  require(future.apply)
  
  # Create output template
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  
  # Apply the test
  if(is.null(data.batch)) { # test without covariats
    data.colData <- factor(data.group)
    testPvals <- future_sapply(X = features,
                               FUN = function(x) {
                                 geneExpr = data.parsed[x, ]
                                 gr1 = geneExpr[data.colData == levels(data.colData)[1]]
                                 gr2 = geneExpr[data.colData == levels(data.colData)[2]]
                                 DifferentialLRT(x = unlist(gr1), y = unlist(gr2))
                               })
    data.out[features, "pval"] <- testPvals
    data.out[features, "FDR"] <- p.adjust(p = data.out[features, "pval"],
                                          method = "fdr")
    if(compute_FC) data.out[, "logFC"] <- total.diff
    data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
    data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
  } else {
    error.json("bimod method doesn't support covariates")
  }
} 

# [Seurat] Run negbinom / poisson by [MASHA] ----------------------------------
# Seurat function DiffExpTest as a base
if(std_method_name == "negbinom_seurat" || 
   std_method_name == "poisson_seurat") {
  if(std_method_name == "negbinom_seurat") std_method_name <- "negbinom-seurat" 
  if(std_method_name == "poisson_seurat") std_method_name <- "poisson-seurat" 
  
  # Check
  if(!is_count_table) {
    error.json(paste0("Data should be a count table in order to run ",std_method_name))
  }
  
  require(Seurat)
  require(future.apply)
  require(MASS)
  
  # info about cell groups
  data.colData <- data.frame(group = factor(data.group))
  # add covariates, if any
  if(!is.null(data.batch)) { # test with covariats
    data.colData <- cbind(data.colData, data.batch)
  }
  
  # apply tests
  testPvals <- future_sapply(X = features,
                             FUN = function(x) applyNegbinomOrPoisson(std_method_name,
                                                                      x, min_pct,
                                                                      data.parsed[x, ], 
                                                                      data.colData))
  testPvals <- unlist(testPvals)
  
  # Create output template
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) <- c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features, "pval"] <- testPvals
  data.out[features, "FDR"] <- p.adjust(p = data.out[features, "pval"],
                                        method = "fdr")
  if(compute_FC) data.out[, "logFC"] <- total.diff
  data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
}

# [Seurat] Run MAST by [MASHA] ------------------------------------------------
if(std_method_name == "mast_seurat") {  
  # Check
  if(!is_log2) {
    error.json("Data should be log in order to run MAST")
  }
  require(Seurat)
  require(MAST)
  
  # info about cell groups
  data.colData <- data.frame(group = factor(data.group))
  # add covariates, if any
  if(!is.null(data.batch)) { # test with covariats
    # scale numeric covariats
    data.batch.scaled <- data.frame(data.batch)
    for (i in 1:ncol(data.batch.scaled)) {
      if (class(data.batch.scaled[, i]) == 'numeric') {
        data.batch.scaled[, i] <- scale(data.batch.scaled[, i])
      } 
    }
    data.colData <- cbind(data.batch.scaled, data.colData)
  }
  
  fdat <- data.frame(primerid = features)
  rownames(fdat) <- features
  
  # perform de with mast
  sca <- MAST::FromMatrix(exprsArray = as.matrix(data.parsed[features, ]),
                          cData = data.colData, fData = fdat)
  # create formula:
  modelAsStr <- paste("~  ", paste0(colnames(data.colData), collapse = ' + '))
  modelEqv <- formula(modelAsStr)
  zlmCond <- MAST::zlm(formula = modelEqv, sca = sca)
  # get p-values
  summaryCond <- summary(object = zlmCond, doLRT = 'group2')
  summaryDt <- summaryCond$datatable
  testPvals <- unlist(summaryDt[summaryDt$component == "H", 4])
  # features order isn't the same as we submitted!
  features.return <- unlist(summaryDt[summaryDt$component == "H", 1])
  features.return <- as.integer(features.return)
  
  # Create output template
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features.return, "pval"] <- testPvals
  data.out[features.return, "FDR"] <- p.adjust(p = data.out[features.return, "pval"],
                                               method = "fdr")
  if(compute_FC) data.out[, "logFC"] <- total.diff
  data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
}

# [Seurat] Run LR by [MASHA] --------------------------------------------------
if(std_method_name == "lr_seurat") {
  require(lmtest)
  require(Seurat)
  require(future.apply)
  
  # info about cell groups
  data.colData <- data.frame(group = factor(data.group))
  # add covariates, if any
  if(!is.null(data.batch)) { # test with covariats
    data.colData <- cbind(data.colData, data.batch)
  }
  
  testPvals <- future_sapply(X = features,
                             FUN = function(x) {
                                  model.data <- cbind(GENE = unlist(data.parsed[x, ]), 
                                                      data.colData)
                                  if (is.null(data.batch)) {
                                     fmla <- as.formula("group ~ GENE")
                                     fmla2 <- as.formula("group ~ 1")
                                  } else {
                                     fmla <- as.formula(paste( "group ~ GENE +",
                                                               paste(colnames(x = data.batch),
                                                                     collapse = "+")))
                                     fmla2 <- as.formula(paste("group ~", 
                                                              paste(colnames(x = data.batch),
                                                                     collapse = "+")))
                                  }
                                  model1 <- glm(formula = fmla, "binomial",
                                                model.data)
                                  model2 <- glm(formula = fmla2, "binomial",
                                                model.data)
                                  lrtest <- lrtest(model1, model2)
                                  return(lrtest$Pr[2])})
  # Create output template
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features, "pval"] <- testPvals
  data.out[features, "FDR"] <- p.adjust(p = data.out[features, "pval"],
                                        method = "fdr")
  if(compute_FC) data.out[, "logFC"] <- total.diff
  data.out[, "AveG1"] <- rowMeans(data.parsed[, cells_1, drop = F])
  data.out[, "AveG2"] <- rowMeans(data.parsed[, cells_2, drop = F])
}

#[Seurat] Run roc by [MASHA] -------------------------------------------------
#Implementation of ROC from seurat, but not functional one, as it doesn't 
#return DE genes
# if(std_method_name == "roc_seurat") {
  # require(ROCR)
  # require(Seurat)
  
#  columns in which different cell groups are residing
  # data.group.fact <- factor(data.group)
  # gr1Cols <- data.group.fact == levels(data.group.fact)[1]
  # gr2Cols <- data.group.fact == levels(data.group.fact)[2]
#  covariates aren't supported
  # if(!is.null(data.batch)) { # test with covariats
    # error.json("Covariates are not supported by roc")
  # }
  
  # to.return <- AUCMarkerTest(data1 = data.parsed[, gr1Cols, drop = F],
                             # data2 = data.parsed[, gr2Cols, drop = F],
                             # mygenes = features)
  # to.return$power <- 2 * abs(to.return$myAUC - 0.5)
  # to.return
# }

# SANITY CHECK for FC direction -----------------------------------------------
not_NA <- !is.na(data.out[,"logFC"])
data.correct <- length(which(sign(data.out[not_NA, "AveG1"] - 
                                  data.out[not_NA, "AveG2"]) == 
                               sign(data.out[not_NA, "logFC"])))
data.incorrect <- sum(not_NA) - data.correct
if(data.correct > data.incorrect) {
  print("SANITY CHECK: OK")
} else {
  # SANITY CHECK (because these DE methods are wild and we want G1 as reference)
  print("SANITY CHECK: ERROR")
  data.out[, "logFC"] <- -data.out[,"logFC"] 
}

# Output results --------------------------------------------------------------
# Open Loom in writing mode for writing results
data.loom <- open_with_lock(input_matrix_filename, "r+")
add_matrix_dataset(handle = data.loom, 
                   dataset_path = output_matrix_dataset,
                   dataset_object = t(data.out))
close_all()

# Volcano plot
if(compute_FC){
  # Prepare text description for each dot:
  not_NA <- !is.na(data.out[,"logFC"])
  data.volcano <- as.data.frame(data.out[not_NA, c("logFC", "pval")])
  data.ens <- data.ens[not_NA]
  data.gene <- data.gene[not_NA]
  data.volcano$text <- paste("log FC: ", round(data.volcano[,"logFC"], 3),
                             "<br>p-value: ", formatC(data.volcano[,"pval"],
                                                      format = "e", digits = 3),
                             "<br>Ensembl: ", data.ens, "<br>Gene: ", data.gene, 
                             sep = "")
  
  p <- ggplot(data.volcano, aes(x = logFC, y = -log10(pval), text = text)) + 
       geom_point(shape = 16) + labs(title = "Volcano plot")

  write(serialize(ggplotly(p, tooltip = "text")), 
        file = paste0(output_dir,"/output.plot.json"), append = F)
}

# Generate default JSON file --------------------------------------------------
stats <- list()
stats$time_idle <- time_idle
# Prepare metadata report
stats$metadata <- list(list(name = output_matrix_dataset, on = "GENE", type = "NUMERIC", nber_cols = 5, nber_rows = nrow(data.out), headers = c("log Fold-Change","p-value","FDR","Avg. Exp. Group 1","Avg. Exp. Group 2")))
if(!is.null(data.warnings)) stats$warnings = data.warnings
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
