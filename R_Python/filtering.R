### ASAP Filtering script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
require(jsonlite)

### Default Parameters
set.seed(42)
input.file <- args[1]
output.folder <- args[2]
algorithm <- args[3]

### Load data
#input.file = "C:/Users/Vincent Gardeux/Dropbox/ASAP/TestData/mASC/ASC.htseq-bow-mmseqEnsg84-runCut-20161121-161121.counts.STAR.txt"
data.parsed = read.table(input.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)
data.warnings = NULL
data.errors = NULL

### Plots that will be generated
ycol <- "Nb Expressed Genes [count > 0]"
expressed.genes.per.sample = data.parsed
expressed.genes.per.sample[expressed.genes.per.sample > 0] = 1
data.plots <- data.frame(name="boxplot.png",description="Gene expression distribution in each sample after filtering.")
data.plots <- rbind(data.plots, data.frame(name="barplot.png", description="Number of genes detected as expressed in each sample."))
data.plots <- rbind(data.plots, data.frame(name="expressed.png", description="Cumulative number of expressed genes when taking cells as bulk [ranked by increasing order]."))

### Functions
error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)
  stop(displayed)
}

### Run Filtering algorithms
if (algorithm == "expressed"){ # default [50]
  percent.expressed = as.double(args[4])
  if(percent.expressed <= 0 || percent.expressed > 100) error.json("Percentage should be between 0 and 100")
  data.sum = rowSums(data.parsed)
  percent.expressed = percent.expressed / 100 # Percent of genes to keep
  threshold = sort(data.sum, decreasing = T)[floor(length(data.sum) * percent.expressed)]
  data.out = data.parsed[data.sum > threshold,]
} else if (algorithm == "coeffofvar"){ # default [50]
  percent.expressed = as.double(args[4])
  if(percent.expressed <= 0 || percent.expressed > 100) error.json("Percentage should be between 0 and 100")
  data.cv <- abs(apply(data.parsed,1,sd) / apply(data.parsed,1,mean))
  data.cv[which(is.nan(data.cv))] = 0
  percent.expressed = percent.expressed / 100 # Percent of genes to keep
  threshold = sort(data.cv, decreasing = T)[floor(length(data.cv) * percent.expressed)]
  data.out = data.parsed[data.cv > threshold,]
} else if (algorithm == "var"){ # default [50]
  percent.expressed = as.double(args[4])
  if(percent.expressed <= 0 || percent.expressed > 100) error.json("Percentage should be between 0 and 100")
  data.var <- apply(data.parsed,1,var)
  data.var[which(is.nan(data.var))] = 0
  percent.expressed = percent.expressed / 100 # Percent of genes to keep
  threshold = sort(data.var, decreasing = T)[floor(length(data.var) * percent.expressed)]
  data.out = data.parsed[data.var > threshold,]
} else if (algorithm == "pagoda"){ # Default [1000,1,1]
  require("scde")
  min.lib.size = as.double(args[4])
  if(min.lib.size < 0) error.json("min.lib.size should be greater than 0")
  min.reads = as.double(args[5])
  if(min.reads < 0) error.json("min.reads should be greater than 0")
  min.detected <- as.double(args[6])
  if(min.detected < 0) error.json("min.detected should be greater than 0")
  if(min.detected > ncol(data.parsed)) error.json(paste0("'Min Detected' should be smaller than the total number of cells/samples i.e. <=",ncol(data.parsed)))
  tryCatch({
    data.out <- as.data.frame(clean.counts(data.parsed, min.lib.size = min.lib.size, min.reads = min.reads, min.detected = min.detected))
  }, error = function(err) {
    print(err)
    error.json("The output matrix is empty. You should put less stringent thresholds. It could be that 'Min library size' is set too large.")
  })
  if(nrow(data.out) == 0) error.json("The output matrix has no more genes. You should put less stringent thresholds. It could be that 'Min Detected' is set too large.")
} else if (algorithm == "scanupc"){ # default [0.95, 1, 1]
  UPC.TRY <- function(data.vector, sample, names){
    tryCatch({
      UPC_RNASeq_Single(data.vector, modelType = "nb", featureNames = names, ignoreZeroes = F)
    }, error = function(err) {
      write(paste0(sample, "\n"), file = "data.warnings.tmp", append = T)
      tryCatch({
        UPC_RNASeq_Single(data.vector, modelType = "nb", featureNames = names, ignoreZeroes = T)
      }, error = function(err) {
        write(paste0(sample, "\n"), file = "data.errors.tmp", append = T)
        rep(0, length(data.vector)) # Return only prob = 0
      })
    })
  }
  upc.probability <- as.double(args[4])
  if(upc.probability < 0 || upc.probability > 1) error.json("UPC.probability should be greater than 0 and less than 1")
  nb.cells.detected <- as.double(args[5])
  if(nb.cells.detected < 0) error.json("nb.cells.detected should be greater or equal to 0")
  if(nb.cells.detected > ncol(data.parsed)) error.json(paste0("'Min Detected' should be smaller than the total number of cells/samples i.e. <=",ncol(data.parsed)))
  nb.cores <- as.double(args[6])
  if(nb.cores < 1) error.json("nb.cores should be greater or equal to 1")
  if(file.exists("data.warnings.tmp")) file.remove("data.warnings.tmp")
  if(file.exists("data.errors.tmp")) file.remove("data.errors.tmp")
  if( nb.cores > 1 ) {
    require(snowfall) # Parallel package
    sfInit( parallel=TRUE, cpus=nb.cores ) # Run on nb.cores cores
    cat( "Running in parallel mode on", sfCpus(), "nodes.\n" )
    sfExport("data.parsed")
    sfExport("UPC.TRY")
    sfLibrary(SCAN.UPC)
    data.upc <- as.data.frame(sfSapply(colnames(data.parsed),function(x) UPC.TRY(data.parsed[,x], x, rownames(data.parsed))))
    sfStop() # stop the parallel package
  } else {
    cat( "Running in sequential mode.\n" )
    require("SCAN.UPC")
    data.upc <- as.data.frame(sapply(colnames(data.parsed),function(x) UPC.TRY(data.parsed[,x], x, rownames(data.parsed))))
  }
  if(file.exists("data.warnings.tmp")){
    data.warnings <- sort(as.character(read.table("data.warnings.tmp")[,1]))
    data.warnings <- paste0(length(data.warnings), " sample(s) were run with option 'ignoreZeroes=T' because they were containing too many zeroes for SCAN.UPC to converge [", paste(data.warnings, collapse = ", "), "]")
    file.remove("data.warnings.tmp")
  }
  if(file.exists("data.errors.tmp")){
    data.errors <- sort(as.character(read.table("data.errors.tmp")[,1]))
    data.errors <- paste0(length(data.errors), " sample(s) failed to converge (even with option 'ignoreZeroes=T') [", paste(data.errors, collapse = ", "), "]")
    file.remove("data.errors.tmp")
  }
  data.upc[data.upc > upc.probability] = 1
  data.upc[data.upc != 1] = 0  
  kept <- rownames(data.upc)[rowSums(data.upc) > nb.cells.detected]
  data.out <- data.parsed[kept,]
  expressed.genes.per.sample = data.upc
  ycol <- paste0("Nb Expressed Genes [UPC probability > ", upc.probability,"]")
} else if (algorithm == "cpm"){ # Default [1,4]
  nb.counts.per.cell <- as.double(args[4])
  if(nb.counts.per.cell < 0) error.json("nb.counts.per.cell should be greater than 0")
  nb.cells.detected <- as.double(args[5])
  if(nb.cells.detected < 0) error.json("nb.cells.detected should be greater than 0")
  if(nb.cells.detected > ncol(data.parsed)) error.json(paste0("'Min Detected' should be smaller than the total number of cells/samples i.e. <=",ncol(data.parsed)))
  
  data.cpm <- as.data.frame(apply(data.parsed, 2, function(x) (x/sum(x))*1000000)) 
  expressed.genes.per.sample = data.cpm
  expressed.genes.per.sample[expressed.genes.per.sample > nb.counts.per.cell] = 1
  expressed.genes.per.sample[expressed.genes.per.sample != 1] = 0
  
  data.cpm <- data.cpm[rowSums(data.cpm > nb.counts.per.cell) >= nb.cells.detected, ]
  data.out <- data.parsed[rownames(data.cpm), ]
  ycol <- paste0("Nb Expressed Genes [CPM > ",nb.counts.per.cell,"]")
} else if (algorithm == "scLVM"){ # Default [1,4]
  require(scLVM)
  require(DESeq2)
  fit.model <- args[4]
  if(is.na(fit.model) || fit.model == "") fit.model = 'log'
  ercc.file <- args[5]
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
    png(paste0(output.folder,"/variable.genes.png"), width=500, height=600, type="cairo")
    data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, method = "fdr", threshold = 0.1, fit_type="counts",sfEndo=estimateSizeFactorsForMatrix(data.parsed), sfERCC=sfERCC)
    dev.off()
    data.plots <- rbind(data.plots, data.frame(name="variable.genes.png", description="scLVM fit of technical noise + variable genes"))
    data.out <- as.data.frame(data.parsed[data.variable.genes, ])
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
    png(paste0(output.folder,"/variable.genes.png"), width=500, height=600, type="cairo")
    data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, fit_type = fit.model, threshold = 0.1)
    dev.off()
    data.plots <- rbind(data.plots, data.frame(name="variable.genes.png", description="scLVM fit of technical noise + variable genes"))
    data.out <- as.data.frame(data.parsed[data.variable.genes, ])
  }
} else error.json("This filtering method is not implemented")

if(nrow(data.out) == 0) error.json("The output matrix has no more genes. You should put less stringent thresholds.")

png(paste0(output.folder,"/boxplot.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
boxplot(data.out, outline=F, las=2, ylab="Raw expression / count")
dev.off()

png(paste0(output.folder,"/barplot.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
barplot(colSums(expressed.genes.per.sample), las=2, ylab=ycol)
dev.off()

data.sort = expressed.genes.per.sample[,names(sort(colSums(expressed.genes.per.sample)))]
data.cumsum = c()
for(i in 1:ncol(data.sort))
{
  data.sort[,1] = data.sort[,i] + data.sort[,1]
  data.cumsum = c(data.cumsum, length(which(data.sort[,1] != 0)))
}
png(paste0(output.folder,"/expressed.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
plot(data.cumsum, type="l", xlab="Nb Cells", ylab = ycol)
dev.off

stats = list()
stats$nber_genes = nrow(data.out)
stats$nber_cells = ncol(data.out)
stats$nber_zeros = length(which(data.out == 0))
stats$list_plots = data.plots
stats$with_negatives = !all(data.out >= 0)
stats$warnings = as.list(c(data.warnings, data.errors))
write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)

data.out.cols = colnames(data.out)
data.out$Genes = rownames(data.out)
write.table(data.out[,c("Genes",data.out.cols)], file=paste0(output.folder,"/output.tab") , sep="\t", row.names = F, quote=F, append = F)