# Integration of the single cell datasets with use of Seurat for ASAP ---------

# Libraries -------------------------------------------------------------------
source("hdf5_lib.R")
require(data.table)
require(jsonlite)
require(loomR)
require(scater)
require(Seurat)
require(patchwork)

# Functions -------------------------------------------------------------------
#' guessSpecie
#' Guesses specie from LOOM file based on submitted specie slot and gene 
#' accession codes
#' @param loomPath string, path to loom file
#' @param specieSlot string, slot name in loom file which should contain specie
#' @param accessSlot string, slot name in loom file which should contain gene
#'                   accession codes
#' @return data frame with fromSpecieSlot and fromFromGene columns
guessSpecie <- function(loomPath, specieSlot = NA, accessSlot = NA) {
  loomData <- open_with_lock(loomPath, "r")
  # try to access specie slot
  specieFromSlot <- NA
  if (!is.na(specieSlot)) {
    specieFromSlot <- unique(fetch_dataset(loomData, specieSlot))
    specieFromSlot <- specieFromSlot[specieFromSlot != '']
    if (is.null(specieFromSlot)) {
      specieFromSlot <- NA
    }
  }
 
  # try to access gene official IDs
  specieFromGene <- NA
  if (!is.na(accessSlot)) {
    geneIDs <- fetch_dataset(loomData, accessSlot)
    specieFromGene <- unique(gsub('[0-9]', '', geneIDs))
    specieFromGene <- gsub('[^[:alnum:] ].*$', '', specieFromGene)
    specieFromGene <- unique(specieFromGene[specieFromGene != ''])
    if (is.null(specieFromGene)) {
      specieFromGene <- NA
    }
  }
  close_file(loomData)
 
  result <- data.frame(fromSpecieSlot = specieFromSlot, 
                       fromFromGene = specieFromGene)
  result
}

#' getMetaData
#' Extracts info about metadata fields from loom file into data table
#' @param loomPath string, path to loom file
#' @return data table with column name as a key
getMetaData <- function(loomPath) {
  # get file structure
  loomData <- open_with_lock(loomPath, "r")
  metaData <- as.data.table(h5ls(loomData))
  metaData <- metaData[group != '/']
  close_file(loomData)
  # add to colnames, except name column, name of the data file
  colnames(metaData) <- paste(colnames(metaData), 
                              gsub('[.].*$', '', basename(loomPath)), 
                              sep = ':')
  colnames(metaData) <- gsub('name:.*', 'name', colnames(metaData)) 
  # set key for future merge
  setkey(metaData, name)
  metaData
}

#' getALlMetaDataFields
#' Lists all fields of metadata of all loom files
#' @param metaList list of data tables containing metadata extrated from loom
#'                 files
#' @return data table with all meta data fields present in at least one
#'         loom file
getALlMetaDataFields <- function(metaList) {
  for (i in 1:length(metaList)) {
    metaList[[i]] <- metaList[[i]][, 1:2]
    colnames(metaList[[i]])[1] <- 'path'
  }
  mergedMeta <- do.call(rbind, metaList)
  mergedMeta <- mergedMeta[!duplicated(mergedMeta), ]
}

#' overlapMetaData
#' Finds common fields between metadata in loom files
#' @param metaList list of data tables containing metadata extrated from loom
#'                 files
#' @return 
overlapMetaData <- function(metaList) {
  # list all fields in all files
  allfields <- lapply(metaList, function(x) x$name)
  allfields <- unique(unlist(allfields))
  
  # determine, which fields are common between the files
  mergedMeta <- metaList[[1]]
  for (i in 2:length(metaList)) {
    mergedMeta <- merge(mergedMeta, metaList[[i]])
    setkey(mergedMeta, name)
  }
  
  # inform, if certain fields were lost
  lostfields <- setdiff(allfields, mergedMeta$name)
  if (length(lostfields) != 0){
    msg <- paste("Only common fields between files will be merged.",
                 "Fields:", paste(lostfields, collapse = ','), 
                 'will not be included.')
    warning(msg)
  }
  
  mergedMeta
}

#' integrateMatrices
#' Integrates datasets with Seurat
#' Based on https://satijalab.org/seurat/v3.1/integration.html
#' @param loomList list of paths to loom files
#' @return seurat object with two assays: 1) RNA - just count matrices pasted
#'         together and 2) integrated - an integrated data set, containing 
#'         nfeatures genes
integrateMatrices <- function(loomList, ...) {
  allData <- list() # all expression matrices are going to be here
  for (loomPath in loomList) {
    # get the expression matrix, skip.validate = T because we can have custom
    # loom
    loomData <- connect(filename = loomPath, mode = "r", skip.validate = T)
    exprMatr <- as.Seurat(loomData)
    # close loom files when done
    loomData$close_all()
    # normalize 
    exprMatr <- NormalizeData(exprMatr, verbose = F)
    exprMatr <- FindVariableFeatures(exprMatr, ...)
    allData[[length(allData) + 1]] <- exprMatr
    exprMatr <- NULL # clean memory
  }
  dataAnchors <- FindIntegrationAnchors(object.list = allData, dims = 1:30)
  dataIntegrated <- IntegrateData(anchorset = dataAnchors, dims = 1:30)
  dataIntegrated
}

# Options ---------------------------------------------------------------------
options(echo = T)
args <- commandArgs(trailingOnly = T)

# Default Parameters  ---------------------------------------------------------
set.seed(42)

input_matrixes_filenames <- args[1] # this actually should be avector
seurat_selection_method <- args[2] # usually "vst", parameter for FindVariableFeatures
seurat_nFeatures <- args[3] # number of features  for FindVariableFeatures
output_dir <- args[4]
output_loom_name <- args[5]

# field indicating specie
specieLoc <- '/row_attrs/genus_species'
# field indicating accession IDs for genes
geneIDsLoc <- '/row_attrs/Accession'

# DEBUG inputs
# Test data downloaded from http://mousebrain.org/loomfiles_level_L1.html
# input_matrixes_filenames <- c('l1_cortex1.loom', 'l1_cortex2.loom',
#                               'l1_cortex3.loom')
# seurat_selection_method <- 'vst'
# seurat_nFeatures <- 2000
# output_dir <- '.'
# output_loom_name <- 'wowowowo'

data.warnings <- NULL
# Step 1: check that data come from the same specie ---------------------------
# merging between species isn't supported!
dataSpecies <- lapply(input_matrixes_filenames, guessSpecie, specieLoc, 
                      geneIDsLoc)
dataSpecies <- do.call(rbind, dataSpecies)
if (length(unique(dataSpecies$fromSpecieSlot)) != 1 &
    length(unique(dataSpecies$fromFromGene)) != 1) {
  msg <- paste("Merging between species isn't supported!",
               "Please submit datasets from one specie!")
  error.json(msg)
}

# Step 2: catalogue meta-data -------------------------------------------------
# list all possible meta-data fields from the loom files to be merged
allMetaData <- lapply(input_matrixes_filenames, getMetaData)
allMetaFields <- getALlMetaDataFields(allMetaData)

# Step 3: merge count data ----------------------------------------------------
# Be careful! This step takes a gazilion of time!
intergrDS <- integrateMatrices(input_matrixes_filenames, 
                               selection.method = seurat_selection_method, 
                               nfeatures = seurat_nFeatures, verbose = F)
# intergDS <- readRDS('abu.Rds')

# Step 4: output loom file ----------------------------------------------------
# Because I don't know how to create a loom file from scratch in R, I'm going
# to copy the smallest of the input files, change it name to the desired output
# name and populate the fields with the values extracted from the integrated
# dataset.
# check input loom sizes
loomSizes <- sapply(input_matrixes_filenames, function(x) file.info(x)$size)
# path to the output loom file
outputLoomPath <- paste0(output_loom_name, '.loom')
# make a copy of the smallest input loom file
file.copy(input_matrixes_filenames[loomSizes == min(loomSizes)], 
          to = outputLoomPath)

# open file for output
outputLoom <- open_with_lock(outputLoomPath, "r+")

# Output matrices
# get integrated matrix and put it into matrix slot
integratedMatrix <- as.matrix(t(intergDS@assays$integrated@data))
add_matrix_dataset(handle = outputLoom, dataset_path = '/matrix',
                   dataset_object = integratedMatrix)
# add corresponding Gene and CellID, just to be sure that they will be in
# metadata
add_array_dataset(handle = outputLoom, 
                  dataset_path = '/row_attrs/Gene',
                  dataset_object = colnames(integratedMatrix))
add_array_dataset(handle = outputLoom, 
                  dataset_path = '/col_attrs/CellID',
                  dataset_object = rownames(integratedMatrix))

# get matrix of integer counts which were just pasted together, and put it to
# the /layers/pasted_matrix slot
pastedMatrix <- as.matrix(t(intergDS@assays$RNA@data))
add_matrix_dataset(handle = outputLoom, 
                   dataset_path = '/layers/pasted_matrix',
                   dataset_object = pastedMatrix)

# Outputting meta-data
# get all integrated meta data
outputMetaData <- intergDS@meta.data
# check, which fields from all avaible ones were not populated by Seurat
notPopulated <- setdiff(allMetaFields$name, colnames(outputMetaData))
# do not touch Gene or CellID, because we already populated them
notPopulated <- notPopulated[!notPopulated %in% c('Gene', 'CellID')]
# remove notPopulated fields from output loom
existingMeta <- as.data.table(h5ls(outputLoom))
toRemove <- existingMeta[name %in% notPopulated]
# this sorting allows to remove from the deepest layers first
toRemove <- toRemove[order(-group)]
for (i in 1:nrow(toRemove)) {
  add_array_dataset(handle = outputLoom, 
                    dataset_path = paste(toRemove$group[i], toRemove$name[i],
                                         sep = '/'),
                    dataset_object = NA)
}
# populate fields which are avaible from outputMetaData
for (i in 1:ncol(outputMetaData)) {
  fieldToPop <- colnames(outputMetaData)[i]
  metaDataToAdd <- outputMetaData[, i]
  names(metaDataToAdd) <- rownames(outputMetaData)
  add_array_dataset(handle = outputLoom, 
                    dataset_path = paste0('/col_attrs/', fieldToPop),
                    dataset_object = metaDataToAdd)
}
close_all()