### ASAP Clustering script
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
input_matrix_dataset <- args[4]
output_matrix_dataset <- args[5]
data.warnings <- NULL
time_idle <- 0

### Functions
ComputeSNN <- function(nn_ranked, prune) {
  .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn_ranked, prune)
}

RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
  .Call('_Seurat_RunModularityClusteringCpp', PACKAGE = 'Seurat', SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename)
}

# Internal helper function to dispatch to various neighbor finding methods
NNHelper <- function(data, query = data, k, method, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  return(
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = 'nn2', args = args)
      },
      "annoy" = {
        error.json("Not implemented yet")
        #args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
        #do.call(what = 'AnnoyNN', args = args)
      },
      error.json("Invalid method. Please choose one of 'rann', 'annoy'")
    )
  )
}

GroupSingletons <- function(ids, SNN) {
  # identify singletons
  singletons <- c()
  singletons <- names(which(table(ids) == 1))
  singletons <- intersect(unique(ids), singletons)
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(unique(ids))
  cluster_names <- setdiff(cluster_names, singletons)
  connectivity <- vector(mode = "numeric", length = length(cluster_names))
  names(connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
      } else {
        connectivity[j] <- mean(subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(names(connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0) message(paste(length(x = singletons), "singletons identified.", length(x = unique(x = ids)), "final clusters."))
  return(ids)
}

RunLeiden <- function( adj_mat, partition.type = c('RBConfigurationVertexPartition', 'ModularityVertexPartition', 'RBERVertexPartition', 'CPMVertexPartition', 'MutableVertexPartition', 'SignificanceVertexPartition', 'SurpriseVertexPartition'), initial.membership = NULL, weights = NULL, node.sizes = NULL, resolution.parameter = 1) {
  if (!py_module_available(module = 'leidenalg')) {
    stop("Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).")
  }
  # import python modules with reticulate
  leidenalg <- import(module = "leidenalg")
  ig <- import(module = "igraph")
  # convert matrix input
  adj_mat <- as.matrix(x = ceiling(x = adj_mat))
  # convert to python numpy.ndarray, then a list
  adj_mat_py <- r_to_py(x = adj_mat)
  adj_mat_py <- adj_mat_py$tolist()
  # convert graph structure to a Python compatible object
  snn_graph <- ig$Graph$Adjacency(adj_mat_py)
  # compute partitions
  part <- switch(
    EXPR = partition.type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      resolution_parameter = resolution.parameter
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial.membership,
      weights = weights
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      node_sizes = node.sizes,
      resolution_parameter = resolution.parameter
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      node_sizes = node.sizes,
      resolution_parameter = resolution.parameter
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial.membership
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial.membership,
      node_sizes = node.sizes,
      resolution_parameter = resolution.parameter
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial.membership,
      weights = weights,
      node_sizes = node.sizes
    ),
    stop("please specify a partition type as a string out of those documented")
  )
  return(part$membership + 1)
}

FindNeighbors <- function(object, nn.method = "rann", annoy.metric = "euclidean", k.param = 20, prune.SNN = 1/15, nn.eps = 0, do.compute.snn=T) {
  if (is.null(dim(object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", call. = FALSE)
    object <- as.matrix(object)
  }
  if (is.null(rownames(object))) error.json("Please provide rownames (cell names) with the input object")
  n.cells <- nrow(object)
  if (n.cells < k.param) {
    data.warnings <<- list(name=paste0("k.param is larger than the number of cells. Set to ", n.cells - 1, "."), description=paste0("Setting k.param to number of cells - 1"))
    k.param <- n.cells - 1
  }
  # find the k-nearest neighbors for each single cell
  message("Computing nearest neighbor graph")
  nn.ranked <- NNHelper( data = object, k = k.param, method = nn.method, searchtype = "standard", eps = nn.eps, metric = annoy.metric)
  nn.ranked <- nn.ranked$nn.idx
  
  # Compute SNN
  if(do.compute.snn){
    message("Computing SNN")
    snn.matrix <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(snn.matrix) <- rownames(object)
    colnames(snn.matrix) <- rownames(object)
    snn.matrix <- as.Graph(snn.matrix)
    
    return(snn.matrix)
  }
  
  # convert nn.ranked into a Graph
  j <- as.numeric(t(nn.ranked))
  i <- ((1:length(j)) - 1) %/% k.param + 1
  nn.matrix <- as(sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(object), nrow(object))), Class = "Graph")
  rownames(nn.matrix) <- rownames(object)
  colnames(nn.matrix) <- rownames(object)
  return(nn.matrix)
}

FindClusters <- function(object, modularity.fxn = 1, resolution = 0.8, algorithm = "louvain", n.start = 10, n.iter = 10, random.seed = 42) {
  if (algorithm == "louvain") {
    ids <- RunModularityClusteringCpp(SNN = object, modularity = modularity.fxn, resolution = resolution, algorithm = 1, nRandomStarts = n.start, nIterations = n.iter, randomSeed = random.seed, printOutput = T, edgefilename = '')
  } else if (algorithm == "louvain_with_multilevel_refinement") {
    ids <- RunModularityClusteringCpp(SNN = object, modularity = modularity.fxn, resolution = resolution, algorithm = 2, nRandomStarts = n.start, nIterations = n.iter, randomSeed = random.seed, printOutput = T, edgefilename = '')
  } else if (algorithm == "SLM") {
    ids <- RunModularityClusteringCpp(SNN = object, modularity = modularity.fxn, resolution = resolution, algorithm = 3, nRandomStarts = n.start, nIterations = n.iter, randomSeed = random.seed, printOutput = T, edgefilename = '')
  } else if (algorithm == "leiden") {
    suppressPackageStartupMessages(require(reticulate))
    ids <- RunLeiden(adj_mat = object, partition.type = "RBConfigurationVertexPartition", initial.membership = NULL, weights = NULL, node.sizes = NULL, resolution.parameter = resolution)
  } else error.json("algorithm not recognised, it should be [louvain, leiden]")
  
  names(ids) <- colnames(object)
  ids <- GroupSingletons(ids = ids, SNN = object)
  
  return(ids)
}

serialize <- function(widget) {
  htmlwidgets:::toJSON2(widget, pretty=TRUE, digits = 3)
}

fviz_silhouette <- function (sil.obj) {
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(plotly))
  if (inherits(sil.obj, c("eclust", "hcut", "pam", "clara", "fanny"))) df <- as.data.frame(sil.obj$silinfo$widths)
  else if (inherits(sil.obj, "silhouette")) df <- as.data.frame(sil.obj[, 1:3])
  else stop("Don't support an oject of class ", class(sil.obj))
  df <- df[order(df$cluster, -df$sil_width), ]
  df$name <- factor(data.cell_names[as.numeric(rownames(df))], levels = data.cell_names[as.numeric(rownames(df))])
  df$cluster <- as.factor(df$cluster)
  #df$sil_width <- round(df$sil_width, 3)
  p <- ggplot(df, aes(x = name, y = sil_width)) + geom_bar(stat = "identity", aes(fill = cluster)) + 
    labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", "\n Average silhouette width: ", round(mean(df$sil_width), 2))) + 
    ggplot2::ylim(c(NA, 1)) + geom_hline(yintercept = mean(df$sil_width), linetype = "dashed", color = "red")
  p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  data.plotly <- plotly::ggplotly(p)
  
  write(serialize(data.plotly), file = paste0(output_dir, "/output.plot.json"), append=F)
}

sihouette.plot <- function(data.clust.p, data.dist.p){
  suppressPackageStartupMessages(require(cluster))
  data.silho <- apply(as.matrix(data.clust.p), 2, function(x) mean(silhouette(x, data.dist.p)[,3]))
  best.k <- names(which(data.silho == max(data.silho)))
  fviz_silhouette(silhouette(data.clust.p[,best.k], data.dist.p)) # Beware any mix between "k" and k. It works for now
  data.clust.p[,best.k]
}

### Open the existing Loom in read-only mode and recuperate the infos (not optimized for OUT-OF-RAM computation)
to_transpose <- T
if(startsWith(input_matrix_dataset, "/col_attrs")) to_transpose <- F
data.loom <- open_with_lock(input_matrix_filename, "r")
data.parsed <- fetch_dataset(data.loom, input_matrix_dataset, transpose = to_transpose) # If run on dimension reduction (like PCA), then do not transpose the matrix
data.cell_names <- fetch_dataset(data.loom, "/col_attrs/CellID") # This makes everything stay lock forever
close_all()
if(is.null(data.parsed)) error.json("This dataset does not exist in the Loom file")

# Replace NaNs by NA
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
data.parsed[is.nan(data.parsed)] <- NA

### Handle clusters
if(std_method_name != "seurat"){
  if(is.na(args[6]) || args[6] == "" || args[6] == "auto" || args[6] == "0") {
    nbclust <- c(2, min(ncol(data.parsed), 11) - 1)
  } else nbclust <- as.numeric(strsplit(args[6], ":")[[1]])
  if(length(nbclust) != 1 && length(nbclust) != 2) error.json("This number of clusters does not make sense")
  if(length(nbclust) == 1 && nbclust > ncol(data.parsed)) error.json("Cannot have more clusters than cells")
  if(length(nbclust) == 2 && nbclust[2] > ncol(data.parsed)) error.json("Max cluster cannot be greater than number of cells")
  if(length(nbclust) == 2 && nbclust[1] < 2) nbclust[1] <- 2
}

### Run clustering algorithm
if (std_method_name == "kmeans"){
  set.seed(42)
  
  # Parameters
  algorithm <- args[7]
  if(is.null(algorithm) | is.na(algorithm)) algorithm <- "Hartigan-Wong"
  
  if(length(nbclust) == 1){
    data.clust <- kmeans(x = t(data.parsed), centers = nbclust, algorithm = algorithm)
    data.out <- data.clust$cluster
  } else {
    data.clust <- sapply(nbclust[1]:nbclust[2], function(i) kmeans(x = t(data.parsed), centers = i, algorithm = algorithm)$cluster)
    colnames(data.clust) = nbclust[1]:nbclust[2]
    data.dist <- dist(t(data.parsed), method = "euclidean")
    data.out <- sihouette.plot(data.clust, data.dist)
  }
} else if (std_method_name == "hclust"){ # Default [euclidean, ward.D2]
  # Parameters
  dist.method <- args[7]
  if(is.null(dist.method) | is.na(dist.method)) dist.method <- "euclidean"
  
  clust.method <- args[8]
  if(is.null(clust.method) | is.na(clust.method)) clust.method <- "ward.D2"
  
  if(dist.method == "pearson" || dist.method == "spearman") { # Distance matrix
    data.dist <- as.dist(1 - cor(data.parsed, method = dist.method))
  } else {
    data.dist <- dist(t(data.parsed), method = dist.method)
  }
  
  # Generating dendrogram
  data.fit <- hclust(data.dist, method=clust.method)
  
  # Continue with generating clustering output
  if(length(nbclust) == 1){
    data.out <- cutree(data.fit, k = nbclust)
  } else if(length(nbclust) == 2){
    data.clust <- sapply(nbclust[1]:nbclust[2],function(i) cutree(data.fit, k = i))
    colnames(data.clust) = nbclust[1]:nbclust[2]
    data.out <- sihouette.plot(data.clust, data.dist)
  }
} else if (std_method_name == "sc3"){ # Default [8, T]
  # Packages
  suppressPackageStartupMessages(require(SC3))
  suppressPackageStartupMessages(require(SingleCellExperiment))
  
  # Parameters
  nb.cores <- args[7]
  if(is.na(nb.cores) || nb.cores == ""){
    print("No nbcores parameter. Running with 8.")
    nb.cores <- 8
  }
  
  # Check if we don't have non-variant cells
  tryCatch({
    if(!all(apply(data.parsed, 2, var, na.rm=TRUE) != 0)) error.json("Cannot rescale a constant/zero column to unit variance. Some Cells have constant/zero variance, SC3 cannot run. Please consider filtering your dataset for empty cells.")
    if(!all(apply(data.parsed, 1, var, na.rm=TRUE) != 0)) error.json("Cannot rescale a constant/zero column to unit variance. Some Genes have constant/zero variance, SC3 cannot run. Please consider filtering your dataset for low/non-expressed genes.")
  }, error = function(err) {
    if(grepl("missing value where TRUE/FALSE needed", err)) error.json("Your dataset contains NA values. It can come from the Normalization/Scaling step, if this was performed on non-filtered data (rows of 0s for e.g.). Please consider filtering your dataset for empty cells and low/non-expressed genes prior normalization/scaling.")
    error.json(err)
  })
  
  # Prepare the data (SCE format) assuming the data is log2, to prevent SC3 to further log it
  data.sc3 <- SingleCellExperiment(assays = list(logcounts = as.matrix(data.parsed)))
  rowData(data.sc3)$feature_symbol = 1:nrow(data.parsed)
  
  # Starting SC3
  best.k <- 2
  if(length(nbclust) == 1){
    data.sc3 <- sc3(data.sc3, ks=nbclust, kmeans_iter_max = 1E4, n_cores=as.integer(nb.cores), gene_filter = F, rand_seed = 42)
    best.k <- as.numeric(nbclust)
  } else if(length(nbclust) == 2){
    data.sc3 <- sc3(data.sc3, ks=nbclust[1]:nbclust[2], kmeans_iter_max = 1E4, n_cores=as.integer(nb.cores), gene_filter=F, rand_seed = 42)
    data.silho <- sapply(nbclust[1]:nbclust[2], function(k) mean(eval(bquote(`$`(data.sc3@metadata$sc3$consensus, .(as.name(k)))))$silhouette[,"sil_width"])) # Get the silhouette width for k in data.sc3@sc3$consensus$`k`$silhouette
    names(data.silho) <- nbclust[1]:nbclust[2]
    best.k <- as.numeric(names(which(data.silho == max(data.silho))))
    fviz_silhouette(data.sc3@metadata$sc3$consensus[[as.character(best.k)]]$silhouette)
  }
  if(ncol(data.parsed) > 5000) data.sc3 <- sc3_run_svm(data.sc3, ks = best.k) # Impute the others
  data.out <- colData(data.sc3)[[paste0("sc3_",best.k,"_clusters")]]
  names(data.out) <- NULL
  #data.out <- eval(bquote(`$`(data.sc3@metadata$sc3$consensus, .(as.name(best.k)))))$silhouette[,"cluster"]
} else if (std_method_name == "seurat"){
  # Packages
  suppressPackageStartupMessages(require(RANN)) # For SNN
  suppressPackageStartupMessages(require(Matrix)) # For sparseMatrix() function
  suppressPackageStartupMessages(require(Seurat))
  
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013)
  
  # Parameters
  k.param <- args[6] # Defines k for the k-nearest neighbor algorithm
  if(is.null(k.param) | is.na(k.param)) k.param <- 20
  else k.param = as.numeric(k.param)
  
  resolution <- args[7] # Use a resolution above (below) 1.0 if you want to obtain a larger (smaller) number of communities/clusters.
  if(is.null(resolution) | is.na(resolution)) resolution <- 0.8
  else resolution = as.numeric(resolution)
  
  algorithm <- args[8] # Algorithm for modularity optimization (louvain; louvain_with_multilevel_refinement; SLM Smart Local Moving; leiden)
  if(is.null(algorithm) | is.na(algorithm)) algorithm <- "louvain"
  
  is.snn <- args[9] # Defines k for the k-nearest neighbor algorithm
  if(is.null(is.snn) | is.na(is.snn)) is.snn <- "snn"
  if(is.snn == "snn"){
    is.snn <- T
  } else if(is.snn == "nn") {
    is.snn <- F
  } else error.json(paste("Last parameter should be 'nn' or 'snn', you entered", is.snn))
    
  # Run first SNN (Shared Nearest Neighbor graph)
  nn.graph <- FindNeighbors(t(data.parsed), k.param = k.param, prune.SNN = 1/15, nn.eps = 0, do.compute.snn = is.snn)
  
  # Run graph-based clustering
  data.out <- FindClusters(nn.graph, resolution = resolution, algorithm = algorithm, n.start = 10, n.iter = 10)
  data.out <- as.numeric(data.out)

  # Correct if 0
  if(min(data.out) == 0) data.out <- data.out + 1
} else error.json("This clustering method is not implemented.")

# Open Loom in writing mode for writing results
data.loom <- open_with_lock(input_matrix_filename, "r+") # I have no clue why it's not properly closed here
add_array_dataset(handle = data.loom, dataset_path = output_matrix_dataset, dataset_object = data.out, storage.mode_param = "integer")
close_all()

# Generate default JSON file
clusts = as.list(table(data.out))
stats <- list()
stats$time_idle <- time_idle
stats$nber_clusters = length(clusts)
# Prepare metadata report
stats$metadata = list(list(name = output_matrix_dataset, on = "CELL", type = "DISCRETE", nber_cols = length(data.out), nber_rows = 1, categories = clusts))
if(!is.null(data.warnings)) stats$warnings = list(data.warnings)
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
