### ASAP Clustering script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Default Parameters
set.seed(42)
input.file <- args[1]
output.folder <- args[2]
algorithm <- args[3]

### Libraries
require(jsonlite)

### Plots that will be generated
data.plots = NULL

### Functions
sihouette.plot <- function(data.clust.p, data.dist.p){
  require(cluster)
  data.silho <- apply(as.matrix(data.clust.p), 2, function(x) mean(silhouette(x, data.dist.p)[,3]))
  best.k <- names(which(data.silho == max(data.silho)))
  data.plots <<- data.frame(name="silhouette.png",description=paste0("Silhouette plot of estimated best k = ",best.k))
  png(paste0(output.folder,"/silhouette.png"), width=500, height=600, type="cairo")
  plot(silhouette(data.clust.p[,best.k], data.dist.p)) # Beware any mix between "k" and k. It works for now
  dev.off()
  data.clust.p[,best.k]
}

error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)
  stop(displayed)
}

### Read file
if(substring(input.file, nchar(input.file) - 4) == ".json") { # if JSON input
  data.json <- fromJSON(paste(readLines(input.file), collapse=""))
  data.norm <- data.frame(text = data.json$text)
  if(!is.null(data.json$PC1)) data.norm$PC1 = data.json$PC1
  if(!is.null(data.json$PC2)) data.norm$PC2 = data.json$PC2
  if(!is.null(data.json$PC3)) data.norm$PC3 = data.json$PC3
  if(!is.null(data.json$PC4)) data.norm$PC4 = data.json$PC4
  if(!is.null(data.json$PC5)) data.norm$PC5 = data.json$PC5
  row.names(data.norm) <- data.norm$text
  data.norm <- data.norm[,-which(colnames(data.norm) %in% c("text"))]
  data.norm <- data.frame(t(data.norm), check.names = F)
} else data.norm <- read.table(input.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)

### Handle clusters
if(is.na(args[4]) || args[4] == "" || args[4] == "auto" || args[4] == "0") {
  nbclust <- c(2, min(ncol(data.norm), 21) - 1)
} else nbclust <- as.numeric(strsplit(args[4], ":")[[1]])
if(length(nbclust) != 1 && length(nbclust) != 2) error.json("This number of clusters does not make sense")
if(length(nbclust) == 1 && nbclust > ncol(data.norm)) error.json("Cannot have more clusters than cells")
if(length(nbclust) == 2 && nbclust[2] > ncol(data.norm)) error.json("Max cluster cannot be greater than number of cells")

### Run clustering algorithm
if (algorithm == "kmeans"){ # default []
  if(length(nbclust) == 1){
    data.clust <- kmeans(t(data.norm), nbclust)
    data.out <- data.clust$cluster
  } else if(length(nbclust) == 2){
    if(nbclust[1] < 2) nbclust[1] <- 2
    data.clust <- sapply(nbclust[1]:nbclust[2],function(i) kmeans(t(data.norm), i)$cluster)
    colnames(data.clust) = nbclust[1]:nbclust[2]
    data.dist <- dist(t(data.norm))
    data.out <- sihouette.plot(data.clust, data.dist)
  }
} else if (algorithm == "pam"){
  require(cluster)
  if(length(nbclust) == 1){
    data.clust <- pam(t(data.norm), nbclust)
    data.out <- data.clust$cluster
  } else if(length(nbclust) == 2){
    if(nbclust[1] < 2) nbclust[1] <- 2
    data.clust <- sapply(nbclust[1]:nbclust[2],function(i) pam(t(data.norm), i)$cluster)
    colnames(data.clust) = nbclust[1]:nbclust[2]
    data.dist <- dist(t(data.norm))
    data.out <- sihouette.plot(data.clust, data.dist)
  }
} else if (algorithm == "hclust"){ # Default [euclidean, ward.D2]
  require(networkD3)
  dist.method <- args[5]
  if(is.na(dist.method) || dist.method == ""){
    print("No dist.method parameter. Running with euclidean.")
    dist.method <- "euclidean"
  }
  clust.method <- args[6]
  if(is.na(clust.method) || clust.method == ""){
    print("No clust.method parameter. Running with ward.D2.")
    clust.method <- "ward.D2"
  }
  if(dist.method == "pearson" || dist.method == "spearman") { # Distance matrix
    data.dist <- as.dist(1 - cor(data.norm, method = dist.method))
  } else {
    data.dist <- dist(t(data.norm), method=dist.method)
  }
  data.fit <- hclust(data.dist, method=clust.method) # Generating dendrogram
  # Create HTML file with Dendrogram
  data.network <- dendroNetwork(data.fit, height = 1000)
  saveNetwork(network = data.network, file = paste0(output.folder,"/output.hclust.html"))
  # Continue with generating clustering output
  if(length(nbclust) == 1){
    data.out <- cutree(data.fit, k = nbclust)
  } else if(length(nbclust) == 2){
    if(nbclust[1] < 2) nbclust[1] <- 2
    data.clust <- sapply(nbclust[1]:nbclust[2],function(i) cutree(data.fit, k = i))
    colnames(data.clust) = nbclust[1]:nbclust[2]
    data.dist <- dist(t(data.norm))
    data.out <- sihouette.plot(data.clust, data.dist)
  }
} else if (algorithm == "sc3"){ # Default [8, T]
  require(SC3)
  require(scater)
  nb.cores <- args[5]
  if(is.na(nb.cores) || nb.cores == ""){
    print("No nbcores parameter. Running with 8.")
    nb.cores <- 8
  }
  is.log2 <- args[6]
  if(is.na(is.log2) || is.log2 == F || is.log2 == "false") { 
    is.log2 = F
    print("The data is assumed to NOT be logged. SC3 will logged it.")
  } else is.log2 = T
  if(!all(data.norm >= 0)) {
    print("Some data are negative. Not logging the data.")
    is.log2 = T
  }
  if(!all(apply(data.norm, 2, var, na.rm=TRUE) != 0)) error.json("Cannot rescale a constant/zero column to unit variance. Some Cells have constant/zero variance, SC3 cannot run.")
  data.sc3 <- newSCESet(exprsData = data.norm)
  data.sc3@logged = is.log2
  best.k <- 2
  if(length(nbclust) == 1){
    data.sc3 <- sc3(data.sc3, ks=nbclust, kmeans_iter_max = 1E4, n_cores=as.integer(nb.cores), gene_filter = F)
    best.k <- as.numeric(nbclust)
  } else if(length(nbclust) == 2){
    if(nbclust[1] < 2) nbclust[1] <- 2
    data.sc3 <- sc3(data.sc3, ks=nbclust[1]:nbclust[2], kmeans_iter_max = 1E4, n_cores=as.integer(nb.cores), gene_filter=F)
    data.silho <- sapply(nbclust[1]:nbclust[2], function(k) mean(eval(bquote(`$`(data.sc3@sc3$consensus, .(as.name(k)))))$silhouette[,"sil_width"])) # Get the silhouette width for k in data.sc3@sc3$consensus$`k`$silhouette
    names(data.silho) <- nbclust[1]:nbclust[2]
    best.k <- as.numeric(names(which(data.silho == max(data.silho))))
    data.plots <- rbind(data.plots, data.frame(name="silhouette.png", description=paste0("Silhouette plot of estimated best k = ",best.k)))
    png(paste0(output.folder,"/silhouette.png"), width=500, height=600, type="cairo")
    sc3_plot_silhouette(data.sc3, k = best.k) # Beware any mix between "k" and k. It works for now
    dev.off()
  }
  data.out <- eval(bquote(`$`(data.sc3@sc3$consensus, .(as.name(best.k)))))$silhouette[,"cluster"]
  names(data.out) <- colnames(data.norm)
} else error.json("This clustering method is not implemented.")

stats <- list()
stats$list_plots = data.plots
stats$clusters = as.list(data.out)
write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output.folder,"/output.json"), append=F)
