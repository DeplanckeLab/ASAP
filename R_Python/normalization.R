### ASAP Normalization script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
require(jsonlite)

### Default Parameters
set.seed(42)
input.file <- args[1]
output.folder <- args[2]
algorithm <- args[3]
batch.file <- args[4]

### Load data
data.parsed <- read.table(input.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)
data.warnings = NULL

### Plots that will be generated
ycol <- ""
data.plots <- data.frame(name="boxplot.norm.png",description="Distribution of the expression of genes in each sample after normalization.")

### Functions
error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)
  stop(displayed)
}

### Run Normalization algorithms
if (algorithm == "scale"){ # default []
  data.out = as.data.frame(scale(data.parsed))
  ycol = "Scaled Expression"  
} else  if (algorithm == "log"){ # default []
  data.out = sign(data.parsed) * log2(1 + abs(data.parsed)) # Take into account normalized data with negative values
  ycol = "Log2 Expression"
} else  if (algorithm == "voom"){ # Default []
  require(limma)
  data.out = as.data.frame(voom(counts=data.parsed, normalize.method="quantile", plot=F)$E)
  ycol = "Log2 Expression (Voom)"
} else if (algorithm == "tmm"){ # default []
  require(edgeR)
  data.dge <- DGEList(counts=data.parsed)
  data.dge <- calcNormFactors(data.dge)
  data.out = as.data.frame(cpm(data.dge, normalized.lib.sizes=TRUE, log = T))
  ycol = "Log2 Expression (TMM / edgeR)"
} else if (algorithm == "rpkm"){ # default [""]
  require(edgeR)
  gene.length.file<- args[5]
  if(is.null(gene.length.file) || is.na(gene.length.file) || gene.length.file == "") error.json("Cannot compute RPKM without gene length information")
  data.gene.length <- read.table(gene.length.file, sep="\t", header=F, row.names=1, check.names=F, stringsAsFactors=F)
  data.dge <- DGEList(counts=data.parsed)
  data.dge <- calcNormFactors(data.dge)
  data.out = as.data.frame(rpkm(data.dge, gene.length = data.gene.length[rownames(data.parsed),1], normalized.lib.sizes=TRUE, log = T))
  ycol = "Expression (RPKM)"
} else if (algorithm == "deseq"){ # default []
  require(DESeq2)
  data.colData = data.frame(row.names=colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design=~1)
  data.dds <- estimateSizeFactors(data.dds)
  data.out <- as.data.frame(counts(data.dds, normalized=TRUE))
  data.out <- sign(data.out) * log2(1 + abs(data.out))
  ycol = "Log2 Expression (DEseq2)"
} else if (algorithm == "vsd"){ # default []
  require(DESeq2)
  data.colData = data.frame(row.names=colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design=~1)
  data.out <- assay(varianceStabilizingTransformation(data.dds, blind=T))
  ycol = "Log2 Expression (varianceStabilizingTransformation - DEseq2)"
} else if (algorithm == "rld"){ # default []
  require(DESeq2)
  data.colData = data.frame(row.names=colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design=~1)
  data.out <- assay(rlogTransformation(data.dds, blind=T))
  ycol = "Log2 Expression (rlogTransformation - DEseq2)"
} else if (algorithm == "scLVM"){ # default []
  require(scLVM)
  require(DESeq2)
  fit.model <- args[5]
  if(is.na(fit.model) || fit.model == "") fit.model = 'log'
  keep.most.variable.genes <- T
  if(is.na(args[6]) || args[6] == "" || args[6] == "null" || args[6] == "false" || args[6] == "0") keep.most.variable.genes = F
  ercc.file <- args[7]
  data.ercc = NULL
  if(!is.null(ercc.file) & !is.na(ercc.file) & ercc.file != "" & ercc.file != "null"){
    data.ercc = read.table(ercc.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)
    if(all(sort(colnames(data.parsed)) %in% sort(colnames(data.ercc)))) {
      print("ERCC file is correct. Process with normalization.")
    } else error.json("ERCC file is NOT correct.")
  }
  if(!is.null(data.ercc)){
    print("Normalize read counts based on ERCC.")
    sfERCC <- NULL
    tryCatch({
      sfERCC <<- estimateSizeFactorsForMatrix(data.ercc[,colnames(data.parsed)])
    }, error = function(err) {
      if(grepl("every gene contains at least one zero", err$message)) error.json("scLVM error: At least one Gene should not contain 0 reads.")
      error.json(err$message)
    })
    data.ercc.sf <- t( t(data.ercc[,colnames(data.parsed)]) / sfERCC )
    data.parsed.sf <- t( t(data.parsed) / sfERCC )
    png(paste0(output.folder,"/tech.noise.fit.png"), width=1000, height=600, type="cairo")
    data.tech.noise <- fitTechnicalNoise(data.parsed.sf,nCountsERCC=data.ercc.sf, fit_type = 'counts')
    dev.off()
    data.plots <- rbind(data.plots, data.frame(name="tech.noise.fit.png", description="scLVM fit of technical noise"))
    if(keep.most.variable.genes){
      png(paste0(output.folder,"/variable.genes.png"), width=500, height=600, type="cairo")
      data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, method = "fdr", threshold = 0.1, fit_type="counts",sfEndo=estimateSizeFactorsForMatrix(data.parsed), sfERCC=sfERCC)
      dev.off()
      data.plots <- rbind(data.plots, data.frame(name="variable.genes.png", description="scLVM fit of technical noise + variable genes"))
      data.out <- as.data.frame(data.parsed.sf[data.variable.genes,])
    } else {
      data.out <- as.data.frame(data.parsed.sf)
    }
    data.out <- log2(1 + data.out)
    ycol = "Log2 ERCC Normalized Counts (scLVM [Brennecke et al.] + DEseq2)"
  } else {
    print("Normalize read counts without ERCC.")
    sfCounts <- NULL
    tryCatch({
      sfCounts <<- estimateSizeFactorsForMatrix(data.parsed)
    }, error = function(err) {
      if(grepl("every gene contains at least one zero", err$message)) error.json("scLVM error: At least one Gene should not contain 0 reads.")
      error.json(err$message)
    })
    data.parsed.sf <- t( t(data.parsed) / sfCounts )
    png(paste0(output.folder,"/tech.noise.fit.png"), width=500, height=600, type="cairo")
    data.tech.noise <- fitTechnicalNoise(data.parsed.sf,use_ERCC = F, fit_type = fit.model)
    dev.off()
    data.plots <- rbind(data.plots, data.frame(name="tech.noise.fit.png", description="scLVM fit of technical noise"))
    if(keep.most.variable.genes){
      png(paste0(output.folder,"/variable.genes.png"), width=500, height=600, type="cairo")
      data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, fit_type = fit.model, threshold = 0.1)
      dev.off()
      data.plots <- rbind(data.plots, data.frame(name="variable.genes.png", description="scLVM fit of technical noise + variable genes"))
      data.out <- as.data.frame(data.parsed.sf[data.variable.genes,])
    } else {
      data.out <- as.data.frame(data.parsed.sf)
    }
    data.out <- log2(1 + data.out)
    if(fit.model == 'log') ycol = "Log2 Normalized Counts (scLVM [log-linear fit] + DEseq2)"
    if(fit.model == 'logvar') ycol = "Log2 Normalized Counts (scLVM [2nd order polynomial regression (loess) fit] + DEseq2)"
  }
} else error.json("This normalization method is not implemented")

png(paste0(output.folder,"/boxplot.norm.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
boxplot(data.out, outline=F, las=2, ylab=ycol)
dev.off()

#BATCH EFFECT CORRECTION
if(!is.null(batch.file) & !is.na(batch.file) & batch.file != "" & batch.file != "null"){
  require(sva)
  print("Batch effect correction requested.")
  print("Reading batch file...")
  data.batch = read.table(batch.file, sep="\t", header=F, row.names=1, check.names=F, stringsAsFactors=F)
  if(all(sort(colnames(data.out)) %in% sort(rownames(data.batch)))) {# If all names are in the batch file (handle filtered cells)
    print("Batch file is correct. Process with ComBat.")
    data.batch = factor(data.batch[colnames(data.out),])
    # Test that no batch has only one element
    if(length(which(table(data.batch) == 1)) != 0) {
      data.warnings <- c(data.warnings, paste0("It's not possible to perform batch effect correction with batches having only one sample [", paste(names(table(data.batch)[table(data.batch) == 1]), collapse = ", "),"]. Batch effect correction was skipped."))
    } else if(length(unique(data.batch)) == 1){
      data.warnings <- c(data.warnings, paste0("It's not possible to perform batch effect correction with only one batch across all samples [", as.character(unique(data.batch)), "]. Batch effect correction was skipped."))
    } else {
        # Test variance in every batch
        for(batch in unique(data.batch)) {
          no.variance.genes = which(apply(data.out[,which(data.batch==batch)], 1, var) == 0)
          if(length(no.variance.genes) != 0){
            data.warnings <- c(data.warnings, paste0(length(no.variance.genes), " gene(s) have no variance (var = 0) in batch ", batch,". ComBat (batch effect correction method) cannot be ran with zero-variance genes. They were removed from the dataset."))
            data.out = data.out[-no.variance.genes,]
          }
        }
        # Run Combat
        png(paste0(output.folder,"/batch.correction.output.png"), width=1000, height=1000, type="cairo")
        data.combat.2 = ComBat(dat=data.out, batch=data.batch, par.prior=T, prior.plots=T)
        dev.off()
        data.plots <- rbind(data.plots, data.frame(name="batch.correction.output.png", description="Output of ComBat for batch effect correction model."))
        data.out = data.combat.2
        png(paste0(output.folder,"/boxplot.norm.after.batch.correction.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
        boxplot(data.out, outline=F, las=2, ylab=paste0(ycol," + batch effect correction [ComBat]"))
        dev.off()
        data.plots <- rbind(data.plots, data.frame(name="boxplot.norm.after.batch.correction.png", description="Distribution of the expression of genes in each sample after normalization & batch effect correction."))
    }
  } else error.json("The batch file is not correct.")
}

stats <- list()
stats$nber_genes = nrow(data.out)
stats$nber_cells = ncol(data.out)
stats$nber_zeros = length(which(data.out == 0))
stats$list_plots = data.plots
stats$warnings = as.list(data.warnings)
write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)

data.out.cols = colnames(data.out)
data.out$Genes = rownames(data.out)
write.table(data.out[,c("Genes",data.out.cols)], file=paste0(output.folder,"/output.tab") , sep="\t", row.names = F, quote=F, append = F)
