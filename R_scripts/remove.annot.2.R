### ASAP DE script
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
require(jsonlite)
require(loomR) # For handling Loom files

### Default Parameters
input_matrix_filename <- args[1]
output_dir <- args[2]
input_matrix_dataset <- args[3]

error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output_dir,"/output.json"), append=F)
  stop(displayed)
}

# Open Loom in writing mode
time_idle = 0
repeat{ # Handle the lock of the file
  isLocked = F
  tryCatch({
    data.loom <- connect(filename = input_matrix_filename, mode = "r+")
    if(data.loom$exists(output_matrix_dataset)) data.loom$link_delete(output_matrix_dataset) # Remove existing dimension reduction with same name
    else data.warnings <<- data.frame(name = "This dataset does not exist in the Loom file.", description = "This dataset does not exist")
    data.loom$close_all()
  }, error = function(err) {
    if(grepl("unable to lock file", err$message)) isLocked = T
    else error.json(err$message)
  })
  if(!isLocked) break
  else {
    message("Sleeping 1sec for file lock....")
    time_idle <<- time_idle + 1
    Sys.sleep(1)
  }
}

# Generate default JSON file
stats <- list()
stats$time_idle <- time_idle
# Prepare metadata report
#stats$metadata = list(list(name = output_matrix_dataset, on = "row", type = "NUMERIC", nber_cols = 5, nber_rows = nrow(data.out)))
if(!is.null(data.warnings)) stats$warnings = as.list(data.warnings)
write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "/output.json"), append=F)
