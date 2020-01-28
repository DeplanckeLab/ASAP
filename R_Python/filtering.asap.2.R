### ASAP Filtering script
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
suppressPackageStartupMessages(require(jsonlite))
suppressPackageStartupMessages(source("hdf5_lib.R"))

### Default Parameters
set.seed(42)
input_matrix_filename <- args[1]
output_dir <- args[2]
std_method_name <- args[3]
toKeep <- NULL
time_idle <- 0

### Functions
serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty=TRUE, digits = 3)
}

rowVars <- function(x) {
  unlist(apply(x, 1, var, na.rm = TRUE))
}

Brennecke <- function (expr_mat, spikes = NA, suppress.plot = FALSE, fdr = 0.1, minBiolDisp = 0.5) {
  fullCountTable <- expr_mat
  if (is.character(spikes)) {
    sp <- rownames(fullCountTable) %in% spikes
    countsSp <- fullCountTable[sp, ]
    countsGenes <- fullCountTable[!sp, ]
  } else if (is.numeric(spikes)) {
    countsSp <- fullCountTable[spikes, ]
    countsGenes <- fullCountTable[-spikes, ]
  } else {
    countsSp <- fullCountTable
    countsGenes <- fullCountTable
  }
  meansSp <- rowMeans(countsSp, na.rm = TRUE)
  varsSp <- rowVars(countsSp)
  cv2Sp <- varsSp/meansSp^2
  meansGenes <- rowMeans(countsGenes, na.rm = TRUE)
  meansGenes[meansGenes <= 0] <- NA
  varsGenes <- rowVars(countsGenes)
  cv2Genes <- varsGenes/meansGenes^2
  minMeanForFit <- unname(quantile(meansSp[which(cv2Sp > 0.3)], 0.8))
  useForFit <- meansSp >= minMeanForFit
  if (sum(useForFit, na.rm = TRUE) < 20) {
    data.warnings[[1]] <- list(name=paste0("Too few spike-ins (",useForFit,") exceed minMeanForFit (=",minMeanForFit,"). Minimum required is 20. Recomputing using all genes."), description=paste0("At least 20 spike-ins are required for fitting."))
    meansAll <- c(meansGenes, meansSp)
    cv2All <- c(cv2Genes, cv2Sp)
    minMeanForFit <- unname(quantile(meansAll[which(cv2All > 0.3)], 0.8))
    useForFit <- meansSp >= minMeanForFit
  }
  if (sum(useForFit, na.rm = TRUE) < 30) {
    data.warnings[[2]] <- list(name=paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."), description=paste0("This message is displayed if less than 30 spike-ins are used for fitting."))
  }
  fit <- statmod::glmgam.fit(cbind(a0 = 1, a1tilde = 1/meansSp[useForFit]), cv2Sp[useForFit])
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])
  psia1theta <- a1
  minBiolDisp <- minBiolDisp^2
  m <- ncol(countsSp)
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- (meansGenes * psia1theta + meansGenes^2 * cv2th)/(1 + cv2th/m)
  p <- 1 - pchisq(varsGenes * (m - 1)/testDenom, m - 1)
  padj <- p.adjust(p, "BH")
  sig <- padj < fdr
  sig[is.na(sig)] <- FALSE
  if (!suppress.plot) {
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(plotly))
    
    if(length(cv2Sp) == length(cv2Genes)) cv2Sp = c()
    if(length(meansSp) == length(meansGenes)) meansSp = c()
    data.plot = data.frame(means=c(meansGenes, meansSp),cv2=c(cv2Genes, cv2Sp),is.ERCC=c(rep(F, length(cv2Genes)), rep(T, length(cv2Sp))))
    if(!is.null(data.ens)) data.plot$ens = c(data.ens, rep("ERCC", length(cv2Sp)))
    if(!is.null(data.gene)) data.plot$gene = c(data.gene, rep("ERCC", length(cv2Sp)))
    # Prepare text description for each dot:
    data.plot$text <- paste0("Means: ", round(data.plot$means, 3), "<br>CV\u00B2: ", round(data.plot$cv2, 3), "<br>Ensembl: ", data.plot$ens, "<br>Gene: ", data.plot$gene, "<br>p-value: ", p, "<br>FDR: ", padj)
    data.plot$colPoint <- c(ifelse(padj < fdr, "black", "darkgrey"), rep("blue", length(cv2Sp)))
    data.plot = subset(data.plot, !is.na(cv2))
    
    p <- ggplot(data.plot, aes(x=means, y=cv2, text=text)) + geom_point(shape = 16, aes(col = colPoint)) + labs(title = "Highly variable genes [Brennecke et al, 2013]")
    p <- p + xlab("Average normalized read count") + ylab("Squared coefficient of variation (CV\u00B2)")
    p <- p + scale_x_log10() + scale_y_log10()
    data.line <- data.frame(xg = 10^seq(log10(min(data.plot$mean, na.rm = T)), log10(max(data.plot$mean, na.rm = T)), length.out = 1000))
    data.line$yg = ((a1)/data.line$xg + a0)
    data.line$text = ""
    p <- p + geom_line(data = data.line, mapping = aes(x = xg, y = yg), col = "#FF000080", lwd = 1)
    
    legend.labels = c("HVG", "Filtered out")
    if(any(data.plot$is.ERCC)) legend.labels = c("HVG", "Spike ins", "Filtered out")
    
    p <- p + scale_colour_identity(name = element_blank(), labels = legend.labels, guide = "legend") + theme(legend.justification = c(1, 1), legend.position = c(1, 1))
          
    write(serialize(ggplotly(p, tooltip = "text")), file = paste0(output_dir,"/output.plot.json"), append=F)
    
    #plot(meansGenes, cv2Genes, xaxt = "n", yaxt = "n", log = "xy", 
    #     xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)", 
    #     col = "white")
    #axis(1, 10^(-2:5), c("0.01", "0.1", "1", "10", "100", "1000", expression(10^4), expression(10^5)))
    #axis(2, 10^(-2:3), c("0.01", "0.1", "1", "10", "100", "1000"), las = 2)
    #abline(h = 10^(-2:1), v = 10^(-1:5), col = "#D0D0D0", lwd = 2)
    #points(meansGenes, cv2Genes, pch = 20, cex = 0.2, col = ifelse(padj < 0.1, "#C0007090", colGenes))
    #if (length(meansSp) < length(meansGenes)) {
    #  points(meansSp, cv2Sp, pch = 20, cex = 0.5, col = colSp)
    #}
    #xg <- 10^seq(-2, 6, length.out = 1000)
    #lines(xg, (a1)/xg + a0, col = "#FF000080", lwd = 3)
    #lines(xg, psia1theta/xg + a0 + minBiolDisp, lty = "dashed", col = "#C0007090", lwd = 3)
  }
  return(sig)
}

DropoutDEM3Drop <- function (fit, ntop = NULL, method = "fdr", qval.thresh = 2, suppress.plot = TRUE) {
  vals <- fit$vals
  coeffs <- NBumiFitDispVsMean(fit, suppress.plot = TRUE)
  exp_size <- exp(coeffs[1] + coeffs[2] * log(vals$tjs/vals$nc))
  droprate_exp <- vector(length = vals$ng)
  droprate_exp_err <- vector(length = vals$ng)
  for (i in 1:vals$ng) {
    mu_is <- vals$tjs[i] * vals$tis/vals$total
    p_is <- (1 + mu_is/exp_size[i])^(-exp_size[i])
    p_var_is <- p_is * (1 - p_is)
    droprate_exp[i] <- sum(p_is)/vals$nc
    droprate_exp_err[i] <- sqrt(sum(p_var_is)/(vals$nc^2))
  }
  droprate_exp[droprate_exp < 1/vals$nc] <- 1/vals$nc
  droprate_obs <- vals$djs/vals$nc
  droprate_obs_err <- sqrt(droprate_obs * (1 - droprate_obs)/vals$nc)
  diff <- droprate_obs - droprate_exp
  combined_err <- sqrt(droprate_exp_err^2 + droprate_obs_err^2)
  Zed <- diff/combined_err
  pvalue <- pnorm(Zed, lower.tail = FALSE)
  names(pvalue) <- names(vals$tjs)
  reorder <- order(pvalue, droprate_exp - droprate_obs)
  out <- pvalue[reorder]
  diff <- diff[reorder]
  qval <- p.adjust(out, method = method)
  if (is.null(ntop)) {
    out <- out[qval < qval.thresh]
    diff <- diff[qval < qval.thresh]
    qval <- qval[qval < qval.thresh]
  } else {
    out <- out[1:ntop]
    diff <- diff[1:ntop]
    qval <- qval[1:ntop]
  }
  outTABLE <- data.frame(Gene = names(out), effect_size = diff, p.value = out, q.value = qval)
  if (!suppress.plot) {
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(plotly))
    data.plot = data.frame(xes = log10(vals$tjs/vals$nc), droprate_obs = droprate_obs)
    data.plot$color = densCols(data.plot$xes, data.plot$droprate_obs, colramp = colorRampPalette(c("grey75", "black")))
    toplot = names(out)
    toplot = names(vals$tjs) %in% toplot
    data.plot$color[toplot] = "darkorange"
    data.plot$droprate_exp = droprate_exp
    
    if(!is.null(data.ens)) data.plot$ens = data.ens[as.numeric(names(fit$sizes))]
    if(!is.null(data.gene)) data.plot$gene = data.gene[as.numeric(names(fit$sizes))]
    # Prepare text description for each dot:
    data.plot$text <- paste("Expression (log10): ", round(data.plot$xes, 3), "<br>Dropout Rate [Observed]: ", round(data.plot$droprate_obs, 3), "<br>Dropout Rate [Expected]: ", round(data.plot$droprate_exp, 3), "<br>Ensembl: ", data.plot$ens, "<br>Gene: ", data.plot$gene, sep="")
    
    data.line = data.frame(xval = data.plot$xes, yval = data.plot$droprate_exp, text = "")
    
    p <- ggplot(data.plot, aes(text = text)) + geom_point(aes(x = xes, y = droprate_obs, col = color)) + theme(legend.position="none") + scale_colour_identity()
    p <- p + geom_line(data = data.line, aes(x = xval, y = yval), col = "dodgerblue" , size = 1) 
    p <- p + xlab("log10(expression)") + ylab("Dropout Rate")
    
    write(serialize(ggplotly(p, tooltip = "text")), file = paste0(output_dir,"/output.plot.json"), append=F)
    
    #xes <- log10(vals$tjs/vals$nc)
    #dens.col <- densCols(xes, droprate_obs, colramp = colorRampPalette(c("grey75", "black")))
    #plot(xes, droprate_obs, col = dens.col, pch = 16, xlab = "", ylab = "")
    #title(ylab = "Dropout Rate", line = 2)
    #title(xlab = "log10(expression)", line = 2)
    #points(xes[toplot], droprate_obs[toplot], col = "darkorange", pch = 16)
    #points(xes, droprate_exp, col = "dodgerblue", pch = 16, cex = 1)
  }
  return(outTABLE)
}

### Load data
# Connect to the loom file in read mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
data.loom <- open_with_lock(input_matrix_filename, "r")
data.parsed <- fetch_dataset(data.loom, "/matrix", transpose = T) # t() because rhdf5 returns the t() of the correct matrix we want
data.ercc <- fetch_dataset(data.loom, "/col_attrs/_ERCCs")
data.ens <- fetch_dataset(data.loom, "/row_attrs/Accession")
data.gene <- fetch_dataset(data.loom, "/row_attrs/Gene")
close_all()

### Run Filtering algorithms
if (std_method_name == "hvg"){
  # Check ERCC
  use_ercc <- args[4]
  if(!is.null(use_ercc) & !is.na(use_ercc) & use_ercc == "true") {
    if(is.null(data.ercc)) error.json("No ERCC found in Loom file") 
    data.parsed <- rbind(data.parsed, data.ercc) # data.ercc should exist
    data.ercc <- (nrow(data.parsed) - nrow(data.ercc) + 1):nrow(data.parsed)
  }
  # Parameters
  fdr_param = args[5]
  if(!is.null(fdr_param) & !is.na(fdr_param)) fdr_param = as.numeric(fdr_param)
  minBiolDisp_param = args[6]
  if(!is.null(minBiolDisp_param) & !is.na(minBiolDisp_param)) minBiolDisp_param = as.numeric(minBiolDisp_param)
  
  # Run HVG
  data.hvg <- Brennecke(expr_mat = data.parsed, spikes = data.ercc, fdr = fdr_param, minBiolDisp=minBiolDisp_param) # spikes = if ERCCs
  
  # If not genes, output error JSON
  if(is.null(data.hvg) || sum(data.hvg) == 0) error.json("The output matrix has no more genes. You should put less stringent thresholds.")
  
  # Preparing JSON for filtering using the Java software
  stats = list() 
  stats$kept_genes = as.numeric(which(data.hvg)) - 1 # These are the indexes to filter
  write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir,"/to_keep.json"), append=F)
  
  # Run Java for filtering the rows in the main matrix and all metadata
  system(paste0("java -jar ASAP.jar -T FilterRows -loom ", input_matrix_filename, " -o ", output_dir, " -m keep -row_names_file ", paste0(output_dir,"/to_keep.json")))
} else if (std_method_name == "m3drop"){
  
  suppressPackageStartupMessages(require(M3Drop))
  
  # Parameters
  fdr_param = args[4]
  if(!is.null(fdr_param) & !is.na(fdr_param)) fdr_param = as.numeric(fdr_param)
  else fdr_param = 0.05
  
  # Filter and generate CPM data
  data.cpm <- NBumiConvertData(data.parsed, is.counts=TRUE)
  
  # Fit DANB model for 10x data
  data.DANB.fit <- NBumiFitModel(data.cpm)
  
  # DE Genes at FDR 5%
  data.DE.Genes <- DropoutDEM3Drop(fit = data.DANB.fit, qval.thresh = fdr_param, suppress.plot=FALSE)
  
  # If not genes, output error JSON
  if(nrow(data.DE.Genes) == 0) error.json("The output matrix has no more genes. You should put less stringent thresholds.")
  
  # Preparing JSON for filtering using the Java software
  stats = list() 
  stats$kept_genes = as.numeric(data.DE.Genes$Gene) - 1
  write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir,"/to_keep.json"), append=F)
  
  # Run Java for filtering the rows in the main matrix and all metadata
  system(paste0("java -jar ASAP.jar -T FilterRows -loom ", input_matrix_filename, " -o ", output_dir, " -m keep -row_names_file ", paste0(output_dir,"/to_keep.json")))
} else error.json("This filtering method is not implemented")

# Check generated JSON file, and add stuff to it
stats <- fromJSON(paste0(output_dir, "/output.json"))
if(!is.null(data.warnings)) stats$warnings = data.warnings # Replace it
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
