### ASAP Heatmap script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Default Parameters
input.file <- args[1]
output.folder <- args[2]
which.dendrograms <- "both"
options(expressions=500000)

### Libraries
require(jsonlite)
require(d3heatmap)
require(pheatmap)

### Functions
error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)
  stop(displayed)
}
serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty=TRUE, digits = 12)
}

### Read file
data.norm <- read.table(input.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)
#max = 500
#data.norm <- read.table("C:/Users/Vincent/Dropbox/ASAP/Scripts/normalization_dl_output.tab", sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)[1:max,]
#output.folder="C:/users/Vincent/Desktop"
data.warnings <- NULL

### Heatmap
if(nrow(data.norm) <= 1 & ncol(data.norm) <= 1){
  error.json(paste0("Input data is too small: [",nrow(data.norm),",",ncol(data.norm),"]"))
}
if(nrow(data.norm) == 1){
  which.dendrograms = "column";
}
if(ncol(data.norm) == 1){
  which.dendrograms = "row";
}

dist.method <- args[3]
if(is.na(dist.method) || dist.method == ""){
  print("No dist.method parameter. Running with euclidean.")
  dist.method <- "euclidean"
}
clust.method <- args[4]
if(is.na(clust.method) || clust.method == ""){
  print("No clust.method parameter. Running with ward.D2.")
  clust.method <- "ward.D2"
}

if(dist.method == "pearson" || dist.method == "spearman") {
  distfun <- function(x) as.dist(1 - cor(x, method = dist.method))
}else{
  distfun <- function(x) dist(x, method=dist.method)
}
hclustfun <- function(x) hclust(x, method=clust.method)

# Heatmap
if(ncol(data.norm) < 1000 && nrow(data.norm) < 1000 && (nrow(data.norm) * ncol(data.norm)) < 100000){
  h = d3heatmap(data.norm, hclustfun = hclustfun, distfun = distfun, dendrogram = which.dendrograms)
  write(serialize(h), file = paste0(output.folder,"/output.heatmap.json"), append=F)
}

# Plotting (I can run the heatmap 4-5 times... I know...)
pdf(paste0(output.folder,"/heatmap.small.pdf"), onefile=FALSE)
pheatmap(data.norm, border_color=NA, hclustfun = hclustfun, distfun = distfun, dendrogram = which.dendrograms)
dev.off()

png(paste0(output.folder,"/heatmap.small.png"))
pheatmap(data.norm, hclustfun = hclustfun, distfun = distfun, dendrogram = which.dendrograms)
dev.off()

pdf(paste0(output.folder,"/heatmap.large.pdf"), height = 5 + nrow(data.norm) / 25, width = 5 + ncol(data.norm) / 25, onefile=FALSE)
pheatmap(data.norm, border_color=NA, fontsize_row = max(12 - 0.4  * sqrt(nrow(data.norm)), 1),  fontsize_col = max(10 - 0.6  * sqrt(ncol(data.norm)), 1), hclustfun = hclustfun, distfun = distfun, dendrogram = which.dendrograms)
dev.off()

png(paste0(output.folder,"/heatmap.large.png"), height = 500 + nrow(data.norm), width = 500 + ncol(data.norm))
pheatmap(data.norm,  fontsize_row = max(12 - 0.4  * sqrt(nrow(data.norm)), 1),  fontsize_col = max(10 - 0.6  * sqrt(ncol(data.norm)), 1), hclustfun = hclustfun, distfun = distfun, dendrogram = which.dendrograms)
dev.off()

stats <- list()
stats$nber_genes = nrow(data.norm)
stats$nber_cells = ncol(data.norm)
write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)

