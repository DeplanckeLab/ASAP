### ASAP DE script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Default Parameters
set.seed(42)
input.file <- args[1]
output.folder <- args[2]
algorithm <- args[3]
batch.file <- args[4]
group.file <- args[5]
type.file <- args[8]

### Libraries
require(jsonlite)

### Functions
error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)
  stop(displayed)
}

### Plots that will be generated
data.plots = NULL

### Load data
data.norm <- read.table(input.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)
data.group <- NULL

### Extra params: cutoffs
p.cutoff <- args[6]
if(is.null(p.cutoff) || is.na(p.cutoff) || p.cutoff == "" || p.cutoff == "null") p.cutoff <- 0.05
fc.cutoff <- args[7]
if(is.null(fc.cutoff) || is.na(fc.cutoff) || fc.cutoff == "" || fc.cutoff == "null") fc.cutoff <- 2

### Batch effect file
data.batch <- NULL
if(!is.null(batch.file) & !is.na(batch.file) & batch.file != "" & batch.file != "null") {
  print("Batch file detected. Will perform batch effect correction.")
  data.batch <- read.table(batch.file, sep="\t", header=F, row.names=1, check.names=F, stringsAsFactors=F)
  if(all(sort(colnames(data.norm)) %in% sort(rownames(data.batch)))) { # If all names are in the batch file (handle filtered cells)
    print("Batch file is correct. Process with it.")
  } else error.json("Batch data is wrong.")
}

### Handle groups
if(is.null(group.file) || is.na(group.file) || group.file == "" || group.file == "null") error.json("Select file is missing")
data.group.json <- fromJSON(group.file)
data.group.1 <- data.group.json$group1
data.group.2 <- data.group.json$group2
if(is.null(data.group.1)) error.json("\"group1\" is not found in JSON select file")
if(!all(data.group.1 %in% colnames(data.norm))) error.json("Samples from Selection 1 were not found in data matrix.")

if(is.null(data.group.2)){
  print("One group file found. Performing this group against all other samples [Marker Genes]")
  data.group <- data.frame(row.names=colnames(data.norm), group = rep(2,ncol(data.norm)))
  data.group[data.group.1,] <- 1
} else {
  print("Two group files found. Performing DE [G1 vs G2] and discarding other samples")
  if(!all(data.group.2 %in% colnames(data.norm))) error.json("Samples from Selection 2 were not found in data matrix.")
  all <- unique(c(data.group.1, data.group.2))
  if(length(all) != length(data.group.1) + length(data.group.2)) error.json("Selections should not overlap")
  data.group <- data.frame(row.names=all, group = rep(1,length(all)))
  data.group[data.group.2,] <- 2
}

if(sum(data.group$group == 1) < 2) error.json("Group 1 should contain at least 2 samples")
if(sum(data.group$group == 2) < 2) error.json("Group 2 should contain at least 2 samples")

# Order by group number
data.group <- data.group[with(data.group, order(group)), , drop=F]
data.norm <- data.norm[,rownames(data.group)]
# Create design matrix from groups
if(is.null(data.batch)){
  data.design <- model.matrix(~data.group[colnames(data.norm),])
} else { # incorporate the batch into the model
  data.design <- model.matrix(~data.group[colnames(data.norm),]+data.batch[colnames(data.norm),])
}

# Run DE algorithms
if(algorithm=="limma"){
  require(limma)
  if(type.file == "original" || type.file == "filtered"){ # if data is a count matrix
    print("Input data is a count table. Voom Normalizing...")
    data.voom <- voom(counts = data.norm, normalize.method="quantile", plot=F)$E # normalize
    fit <- NULL
    tryCatch({
      fit <<- lmFit(data.voom, data.design)
    }, warning = function(war) {
      if(grepl("Partial NA coefficients", war$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
    })
    fit <- eBayes(fit)
    DE <- topTable(fit,n=Inf,adjust="fdr") # Create the table with all DE results
    DE <- DE[which(DE$adj.P.Val <= as.double(p.cutoff) & abs(DE$logFC) >= log2(as.double(fc.cutoff))),] # restrict the table to the 2 given cutoff
    DE <- DE[with(DE, order(abs(logFC), decreasing=T)), ]  ## ON ORDONNE par log fold change, n�gatif = upr�gul� dans le groupe 1
    DE$AvCountSelected <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
    DE$AvCountOthers <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
    DE$AvVoomSelected <- rowMeans(data.voom[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
    DE$AvVoomOthers <- rowMeans(data.voom[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
    data.out <- DE[,c(-2, -3, -6)]
    colnames(data.out) <- c("logFC", "pval", "FDR", "AvCountG1", "AvCountG2", "AvNormG1", "AvNormG2")
  } else { 
    fit <- NULL
    tryCatch({
      fit <<- lmFit(data.norm, data.design)
    }, warning = function(war) {
      if(grepl("Partial NA coefficients", war$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
    })
    fit <- eBayes(fit)
    DE <- topTable(fit,n=Inf,adjust="fdr") # Create the table with all DE results
    DE <- DE[which(DE$adj.P.Val <= as.double(p.cutoff) & abs(DE$logFC) >= log2(as.double(fc.cutoff))),] # restrict the table to the 2 given cutoff
    DE <- DE[with(DE, order(abs(logFC), decreasing=T)), ]  ## ON ORDONNE par log fold change
    DE$AvExprSelected <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
    DE$AvExprOthers <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
    data.out <- DE[,c(-2, -3, -6)]
    colnames(data.out) <- c("logFC", "pval", "FDR", "AvNormG1", "AvNormG2")
  }
} else if(algorithm=="deseq2"){ # Assuming data is a count table
  require(DESeq2)
  require(BiocParallel)
  if(is.null(args[9]) || is.na(args[9]) || args[9] == "" || args[9] == "null") { 
    nb.cores <- 8
  } else {
    nb.cores <- as.numeric(args[9])
  }
  if(is.null(data.batch)){
    data.colData = data.frame(row.names=colnames(data.norm), group=factor(data.group[colnames(data.norm),]))
    data.dds <- DESeqDataSetFromMatrix(data.norm, colData=data.colData, design=~group)
  } else { # incorporate the batch into the model
    data.colData = data.frame(row.names=colnames(data.norm), group=factor(data.group[colnames(data.norm),]), batch=factor(data.batch[colnames(data.norm),]))
    data.dds <- NULL
    tryCatch({
      data.dds <<- DESeqDataSetFromMatrix(data.norm, colData=data.colData, design=~batch+group)
    }, error = function(err) {
      if(grepl("the model matrix is not full rank", err$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
      error.json(err$message)
    })
  }
  if(nb.cores > 1) {
    register(MulticoreParam(nb.cores))
    data.dds <- DESeq(data.dds, parallel=T)
  } else data.dds <- DESeq(data.dds, parallel=F)
  DE <- results(data.dds)
  DE <- DE[which(DE$padj <= as.double(p.cutoff) & abs(DE$log2FoldChange) >= log2(as.double(fc.cutoff))),]
  DE <- data.frame(DE[with(DE, order(abs(log2FoldChange), decreasing = T)), ])
  data.deseq <- as.data.frame(counts(data.dds, normalized=TRUE))
  data.deseq <- sign(data.deseq) * log2(1 + abs(data.deseq))
  DE$AvCountSelected <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
  DE$AvCountOthers <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
  DE$AvDESeqSelected <- rowMeans(data.deseq[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
  DE$AvDESeqOthers <- rowMeans(data.deseq[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
  data.out <- DE[,c(-1, -3, -4)]
  colnames(data.out) <- c("logFC", "pval", "FDR", "AvCountG1", "AvCountG2", "AvNormG1", "AvNormG2")
} else if(algorithm=="edger"){ # Assuming data is a count table
  require(edgeR)
  data.dge <- DGEList(counts=data.norm, group=data.group[colnames(data.norm),1])
  #data.dge <- relevel(data.dge, ref="1")
  data.dge <- calcNormFactors(data.dge)
  if(is.null(data.batch)){
    data.dge <- estimateDisp(data.dge)
  } else {
    tryCatch({
      data.dge <<- estimateDisp(data.dge, design=data.design, robust=TRUE)
    }, error = function(err) {
      if(grepl("Design matrix not of full rank", err$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
      error.json(err$message)
    })
  }
  DE <- exactTest(data.dge)
  DE <- data.frame(topTags(DE, n=Inf, adjust="fdr"))
  DE <- DE[which(DE$FDR <= as.double(p.cutoff) & abs(DE$logFC) >= log2(as.double(fc.cutoff))),]
  DE <- DE[with(DE, order(abs(logFC), decreasing = T)), ]
  data.edger <- as.data.frame(cpm(data.dge, normalized.lib.sizes=TRUE, log = T))
  DE$AvCountG1 <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
  DE$AvCountG2 <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
  DE$AvEdgeRG1 <- rowMeans(data.edger[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
  DE$AvEdgeRG2 <- rowMeans(data.edger[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
  data.out <- DE[,c(-2)]
  colnames(data.out) <- c("logFC", "pval", "FDR", "AvCountG1", "AvCountG2", "AvNormG1", "AvNormG2")
} else if(algorithm=="scde"){ # Assuming data is a count table
  require(scde)
  if(is.null(args[9]) || is.na(args[9]) || args[9] == "" || args[9] == "null") { 
    nb.cores <- 8
  } else {
    nb.cores <- as.numeric(args[9])
  }
  if(!is.null(data.batch)) { # sanity check before the long computing
    require(limma)
    data.voom <- voom(counts = data.norm, normalize.method="quantile", plot=F)$E # normalize
    tryCatch({
      fit <- lmFit(data.voom, data.design)
    }, warning = function(war) {
      if(grepl("Partial NA coefficients", war$message)) error.json("Seems that your batches and selections are confounded. There is no way to properly analyze the experiment in that case.")
    })
  }
  nb.cores <- as.integer(nb.cores)
  data.model.ifm <- scde.error.models(counts=data.norm, min.nonfailed = 1, min.size.entries = 1, groups=data.group[colnames(data.norm),1], n.cores=nb.cores, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
  data.model.ifm <- data.model.ifm[data.model.ifm$corr.a > 0, ] # keep only valid cells
  data.model.prior <- scde.expression.prior(models=data.model.ifm, counts=data.norm[,rownames(data.model.ifm)], length.out=400, show.plot=F)
  if(is.null(data.batch)){
    DE <- scde.expression.difference(data.model.ifm, data.norm[,rownames(data.model.ifm)], data.model.prior, groups=factor(data.group[rownames(data.model.ifm),1]), n.randomizations=100, n.cores=nb.cores, verbose=1)
  } else {
    DE <- scde.expression.difference(data.model.ifm, data.norm[,rownames(data.model.ifm)], data.model.prior, batch=factor(data.batch[rownames(data.model.ifm),1]), groups=factor(data.group[rownames(data.model.ifm),1]), n.randomizations=100, n.cores=nb.cores, verbose=1)
  }  
  DE$pval <- 2*pnorm(-abs(DE$Z)) # Calculation of p-pvalue from Z-score
  DE$padj <- p.adjust(DE$pval, method="fdr") # Multiple test adjustment
  DE <- DE[which(DE$padj <= as.double(p.cutoff) & abs(DE$mle) >= log2(as.double(fc.cutoff))),]
  DE <- DE[with(DE, order(abs(mle), decreasing = T)), ]
  DE$AvCountSelected <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 1],drop=F])
  DE$AvCountOthers <- rowMeans(data.norm[rownames(DE),rownames(data.group)[data.group == 2],drop=F])
  data.out <- DE[, c(-1, -3, -4, -5, -6)]
  colnames(data.out) <- c("logFC", "pval", "FDR", "AvCountG1", "AvCountG2")
}

# TODO Volcano plot

data.neg <- data.out[data.out$logFC < 0,]
l.neg <- length(which(data.neg$AvCountG1 < data.neg$AvCountG2))
l.pos <- length(which(data.neg$AvCountG1 > data.neg$AvCountG2))
if(l.neg == 0 && l.pos == 0) {
  l.neg <- length(which(data.neg$AvNormG1 < data.neg$AvNormG2))
  l.pos <- length(which(data.neg$AvNormG1 > data.neg$AvNormG2))
}

if(l.neg < l.pos) data.out$logFC <- -data.out$logFC # SANITY CHECK (because these DE methods are wild and we want G1 as reference)

data.out.up <- data.out[data.out$logFC >= 0,]
data.out.up$text = rownames(data.out.up)

data.out.down <- data.out[data.out$logFC < 0,]
data.out.down$text = rownames(data.out.down)

write(toJSON(as.list(data.out.up), method="C", auto_unbox=F), file = paste0(output.folder,"/output.up.json"), append=F)
write(toJSON(as.list(data.out.down), method="C", auto_unbox=F), file = paste0(output.folder,"/output.down.json"), append=F)

