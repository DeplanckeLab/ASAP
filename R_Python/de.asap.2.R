### ASAP DE script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
require(jsonlite)
require(loomR) # For handling Loom files
require(ggplot2)
require(plotly)

### Functions
FindMarkers.default <- function(object, cells.1 = NULL, cells.2 = NULL, fc.threshold = 1.3, test.use = "wilcox", min.pct = 0.1, min.diff.pct = -Inf, max.cells.per.ident = Inf, latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3) {
  if (!(test.use %in% c('negbinom', 'poisson', 'MAST', "LR")) && !is.null(x = latent.vars)) error.json("'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests")
  
  # Perform DE
  de.results <- switch(EXPR = test.use,
    'wilcox' = WilcoxDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2),
    'bimod' = DiffExpTest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2),
    'roc' = MarkerTest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2),
    't' = DiffTTest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2),
    'negbinom' = GLMDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2, min.cells = min.cells.feature, latent.vars = latent.vars, test.use = test.use),
    'poisson' = GLMDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2, min.cells = min.cells.feature, latent.vars = latent.vars, test.use = test.use),
    'MAST' = MASTDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2, latent.vars = latent.vars),
    "LR" = LRDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], cells.1 = cells.1, cells.2 = cells.2, latent.vars = latent.vars),
    error.json("Unknown test: ", test.use))
  
  # Shape DE Results the way we need it
  de.results[, "avg_logFC"] <- total.diff[rownames(x = de.results)]
  de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  de.results$p_val_adj = p.adjust(p = de.results$p_val, method = "bonferroni", n = nrow(object)) # TODO correction
  
  if (test.use == "roc") de.results <- de.results[order(-de.results$power, -de.results$avg_logFC), ]
  else de.results <- de.results[order(de.results$p_val, -de.results$avg_logFC), ]

  return(de.results)
}

serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty=TRUE, digits = 3)
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

### Default Parameters
set.seed(42)
input_matrix_filename <- args[1]
output_dir <- args[2]
std_method_name <- args[3]
input_matrix_dataset <- args[4]
output_matrix_dataset <- args[5]
batch_dataset <- args[6]
group_dataset <- args[7]
group_1 <- args[8]
group_2 <- args[9]
if(group_1 == group_2) error.json("Cannot compute DE from the same group.")
is_count_table <- args[10]
if(is_count_table == "true") { 
  is_count_table <- T
} else if(is_count_table == "false") {
  is_count_table <- F
} else {
  error.json("is_count_table should be 'true' or 'false'")
}
min_pct <- as.numeric(args[11])
min_diff_pct <- as.numeric(args[12])
if(is.null(min_diff_pct) || is.na(min_diff_pct) || min_diff_pct == "" || min_diff_pct == "null"){ 
  min_diff_pct <- -Inf
} else {
  min_diff_pct <- as.numeric(min_diff_pct)
}
fc_threshold <- as.numeric(args[13])
max_cells_per_ident <- args[14]
if(is.null(max_cells_per_ident) || is.na(max_cells_per_ident) || max_cells_per_ident == "" || max_cells_per_ident == "null"){ 
  max_cells_per_ident <- Inf
} else {
  max_cells_per_ident <- as.numeric(max_cells_per_ident)
}
data.warnings <- NULL
time_idle <- 0

### Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
data.loom <- open_with_lock(input_matrix_filename, "r")
if(!data.loom$exists(input_matrix_dataset)) error.json("This dataset does not exist in the Loom file")
data.parsed <- t(data.loom[[input_matrix_dataset]][, ]) # t() because loomR returns the t() of the correct matrix we want
if(!data.loom$exists(group_dataset)) error.json(paste0(group_dataset, " does not exist in the Loom file"))
data.groups <- data.loom[[group_dataset]][]
if(data.loom$exists("/row_attrs/Accession")) data.ens <- data.loom[["/row_attrs/Accession"]][]
if(data.loom$exists("/row_attrs/Gene")) data.gene <- data.loom[["/row_attrs/Gene"]][]
if(is.null(batch_dataset) || is.na(batch_dataset) || batch_dataset == "" || batch_dataset == "null"){ 
  data.batch <- NULL
} else {
  message("Covariates detected. Will add covariates to the model.")
  ### TODO: HANDLE COVARIATES [,,,,,]
  #if (!(std_method_name %in% c('negbinom', 'poisson', 'MAST', "LR")) && !is.null(x = latent.vars)) error.json("'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests")
  data.batch <- data.loom[[batch_dataset]][]
}
data.loom$close_all()

### Assessing log2ity of the data (could also be ln or log10 that should still be ok)
is_log2 <- T
compute_FC <- T
if(is_count_table | max(data.parsed) > 30) is_log2 <- F # Arbitrary cutoff, could probably be set lower than that
if(min(data.parsed) < 0) {
  compute_FC <- F
  is_log2 <- F
  is_count_table <- F
  data.warnings[[1]] <- list(name = "Data has negative values, no FC will be computed", description = "Select another dataset if you want to compute FC.")
}

### Handle groups
groups <- unique(data.groups)
data.group <- rep(NA, ncol(data.parsed))
data.group[data.groups == group_1] <- 1
if(group_2 != "null"){
  message("Two group files found. Performing DE [G1 vs G2] and discarding other samples")
  data.group[data.groups == group_2] <- 2
} else {
  print("One group file found. Performing this group against all other samples [Marker Genes]")
  data.group[data.groups != group_1] <- 2
  group_2 <- 2
}
# Filter NA groups
toKeep = which(!is.na(data.group))
data.parsed <- data.parsed[,toKeep]
data.group <- data.group[toKeep]
data.batch <- data.batch[toKeep]
cells_1 <- data.group == 1
cells_2 <- data.group == 2
features = 1:nrow(data.parsed)

### Check if computable
if(sum(data.groups == group_1) < 3) error.json("Group 1 should contain at least 3 samples")
if(sum(data.groups == group_2) < 3) error.json("Group 2 should contain at least 3 samples")

# [Seurat] Feature selection (based on percentages) to speed up computation
if(is_count_table | is_log2){ # Because checking genes > 0 as "detected"
  pct.1 <- round(x = rowSums(x = data.parsed[, cells_1, drop = FALSE] > 0) / sum(x = cells_1), digits = 3)
  pct.2 <- round(x = rowSums(x = data.parsed[, cells_2, drop = FALSE] > 0) / sum(x = cells_2), digits = 3)
  data.alpha <- cbind(pct.1, pct.2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  features <- which(x = alpha.min > min_pct)
  if(length(x = features) == 0) error.json("No features pass min_pct threshold")
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
  features <- which(x = alpha.min > min_pct & alpha.diff > min_diff_pct)
  if(length(x = features) == 0) error.json("No features pass min_diff_pct threshold")
} else data.warnings[[2]] <- list(name = "Data is neither a count matrix nor detected as logged. No pre-filtering (using &#8217;Min Pct&#8217; and &#8217;Min Diff Pct&#8217;) was performed", description = "Select another dataset if you want to perform pre-filtering (because looking for &#8217;detected genes&#8217; i.e. with values of expression &gt; 0).")

# [Seurat] gene selection (based on average difference/FC) to speed up computation
total.diff <- NULL
if(compute_FC){
  if(!is_log2){ # If also !is_count_table then I don't know the data type, so I still log
    data.1 <- apply(X = data.parsed[, cells_1, drop = FALSE], MARGIN = 1, FUN = function(x) {m = mean(x); return(sign(m) * log2(abs(m) + 1))})
    data.2 <- apply(X = data.parsed[, cells_2, drop = FALSE], MARGIN = 1, FUN = function(x) {m = mean(x); return(sign(m) * log2(abs(m) + 1))})
  } else { # should contain negative values
    data.1 <- apply(X = data.parsed[, cells_1, drop = FALSE], MARGIN = 1, FUN = function(x) return(log2(mean(2^x - 1) + 1)))
    data.2 <- apply(X = data.parsed[, cells_2, drop = FALSE], MARGIN = 1, FUN = function(x) return(log2(mean(2^x - 1) + 1)))
  } 
  total.diff <- data.1 - data.2
  features.diff <- which(abs(total.diff) > log2(fc_threshold))
  features <- intersect(x = features, y = features.diff)
  if(length(x = features) == 0) error.json("No features pass 'FC threshold'")
}

# [Seurat] Downsampling cells if too many per group
if (max_cells_per_ident < Inf) {
  set.seed(42)
  if (sum(cells_1) > max_cells_per_ident) cells_1_kept <- sample(x = which(cells_1), size = max_cells_per_ident)
  else cells_1_kept = which(cells_1)
  if (sum(cells_2) > max_cells_per_ident) cells_2_kept <- sample(x = which(cells_2), size = max_cells_per_ident)
  else cells_2_kept = which(cells_2)
  
  # Filter the datasets
  toKeep = sort(c(cells_1_kept, cells_2_kept))
  data.parsed <- data.parsed[, toKeep] # We keep all features
  data.group <- data.group[toKeep]
  data.batch <- data.batch[toKeep]
  cells_1 <- data.group == 1
  cells_2 <- data.group == 2
  #  if (!is.null(x = latent.vars)) latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
}

# Run DE algorithms
if(std_method_name=="limma"){
  # Package
  require(limma)
  
  # Create design matrix from groups
  if(is.null(data.batch)){
    data.design <- model.matrix(~data.group)
  } else { # incorporate the batch into the model
    data.design <- model.matrix(~data.group+data.batch)
  }
  
  # Test if the data is count or normalized already
  if(is_count_table){ # if data is a count matrix
    message("Input data is a count table. Voom Normalizing...")
    data.voom <- voom(counts = data.parsed[features,], normalize.method="quantile", plot=F)$E # normalize
    fit <- NULL
    tryCatch({
      fit <- lmFit(data.voom, data.design)
    }, warning = function(war) {
      if(grepl("Partial NA coefficients", war$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
    })
  } else { 
    fit <- NULL
    tryCatch({
      if(!is_log2 | !compute_FC) fit <- lmFit(data.parsed[features,], data.design)
      else fit <- lmFit(log2(1 + data.parsed[features,]), data.design)
    }, warning = function(war) {
      if(grepl("Partial NA coefficients", war$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
    })
  }
  
  # Run the DE model
  fit <- eBayes(fit)
  DE <- topTable(fit, n=Inf, adjust="fdr", sort.by = "none") # Create the table with all DE results
  
  # Create output
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features, "logFC"] = DE$logFC
  data.out[features, "pval"] = DE$P.Value
  data.out[features, "FDR"] = DE$adj.P.Val
  data.out[, "AveG1"] = rowMeans(data.parsed[, cells_1, drop=F])
  data.out[, "AveG2"] = rowMeans(data.parsed[, cells_2, drop=F])
} else if(std_method_name == "deseq2"){ # Assuming data is a count table
  # Packages
  require(DESeq2)
  require(BiocParallel)
  
  # Check
  if(!is_count_table) error.json("Data should be a count table in order to run DESeq2")
  
  # Parameters
  nb.cores <- args[15]
  if(!is.null(nb.cores) & !is.na(nb.cores) & nb.cores != "" & nb.cores != "null") nb.cores <- as.numeric(nb.cores)
  else nb.cores <- 8

  data.dds <- NULL
  if(is.null(data.batch)){
    data.colData = data.frame(group=factor(data.group, levels = c(1,2)))
    data.dds <- DESeqDataSetFromMatrix(data.parsed[features,], colData=data.colData, design=~group)
  } else { # incorporate the batch into the model
    data.colData = data.frame(group=factor(data.group, levels = c(1,2)), batch=factor(data.batch))
    tryCatch({
      data.dds <- DESeqDataSetFromMatrix(data.parsed[features,], colData=data.colData, design=~batch+group)
    }, error = function(err) {
      if(grepl("the model matrix is not full rank", err$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
      error.json(err$message)
    })
  }
  tryCatch({
    if(nb.cores > 1) {
      register(MulticoreParam(nb.cores))
      data.dds <- DESeq(data.dds, parallel=T)
    } else data.dds <- DESeq(data.dds, parallel=F)
  }, error = function(err) {
    if(grepl("every gene contains at least one zero", err$message)) error.json("Every gene contains at least one zero, cannot compute log geometric means for estimating size factors.")
    if(!is.null(err)) error.json(err$message)
    #error.json(err)
  })
  
  # Obtain the DE results
  DE <- results(data.dds)
  
  # Create output
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features, "logFC"] = DE$log2FoldChange
  data.out[features, "pval"] = DE$pvalue
  data.out[features, "FDR"] = DE$padj
  data.out[, "AveG1"] = rowMeans(data.parsed[, cells_1, drop=F])
  data.out[, "AveG2"] = rowMeans(data.parsed[, cells_2, drop=F])
} else if(std_method_name == "wilcox-seurat"){ # Seurat function
  require(Seurat)
  require(future.apply)
  
  # Design matrix
  data.colData = data.frame(group=factor(data.group, levels = c(1,2)))
  
  # Create output
  data.out <- matrix(nrow = nrow(data.parsed), ncol = 5)
  colnames(data.out) = c("logFC", "pval", "FDR", "AveG1", "AveG2")
  data.out[features, "pval"] = future_sapply(X = features, FUN = function(x) return(wilcox.test(data.parsed[x, ] ~ data.colData[, "group"])$p.value))
  data.out[features, "FDR"] = p.adjust(p = data.out[features, "pval"], method = "fdr")
  if(compute_FC) data.out[, "logFC"] <- total.diff
  data.out[, "AveG1"] = rowMeans(data.parsed[, cells_1, drop=F])
  data.out[, "AveG2"] = rowMeans(data.parsed[, cells_2, drop=F])
} else error.json(paste0(std_method_name, " method is not implemented"))

# SANITY CHECK for FC direction
not_NA <- !is.na(data.out[,"logFC"])
data.correct <- length(which(sign(data.out[not_NA,"AveG1"] - data.out[not_NA,"AveG2"]) == sign(data.out[not_NA,"logFC"])))
data.incorrect <- sum(not_NA) - data.correct
if(data.correct > data.incorrect) {
  print("SANITY CHECK: OK")
} else {
  print("SANITY CHECK: ERROR")
  data.out[,"logFC"] <- -data.out[,"logFC"] # SANITY CHECK (because these DE methods are wild and we want G1 as reference)
}

# Open Loom in writing mode for writing results
data.loom <- open_with_lock(input_matrix_filename, "r+")
if(data.loom$exists(output_matrix_dataset)) data.loom$link_delete(output_matrix_dataset) # Remove existing dimension reduction with same name
data.loom[[output_matrix_dataset]] = t(data.out)
data.loom$close_all()

# Volcano plot
if(compute_FC){
  # Prepare text description for each dot:
  not_NA <- !is.na(data.out[,"logFC"])
  data.volcano <- as.data.frame(data.out[not_NA, c("logFC", "pval")])
  data.ens <- data.ens[not_NA]
  data.gene <- data.gene[not_NA]
  data.volcano$text <- paste("log FC: ", round(data.volcano[,"logFC"], 3), "<br>p-value: ", formatC(data.volcano[,"pval"], format = "e", digits = 3), "<br>Ensembl: ", data.ens, "<br>Gene: ", data.gene, sep="")
  
  p <- ggplot(data.volcano, aes(x = logFC, y = -log10(pval), text = text)) + geom_point(shape = 16) + labs(title = "Volcano plot")

  write(serialize(ggplotly(p, tooltip = "text")), file = paste0(output_dir,"/output.plot.json"), append=F)
}

# Generate default JSON file
stats <- list()
stats$time_idle <- time_idle
# Prepare metadata report
stats$metadata = list(list(name = output_matrix_dataset, on = "GENE", type = "NUMERIC", nber_cols = 5, nber_rows = nrow(data.out)))
if(!is.null(data.warnings)) stats$warnings = data.warnings
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
