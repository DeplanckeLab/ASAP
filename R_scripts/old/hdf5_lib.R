require(rhdf5) # For handling Loom files # loomR is deprecated

error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(output_dir,"/output.json"), append=F)
  close_all()
  stop(displayed)
}

# Close All handles
close_all <- function() {
  h5closeAll()
}

# Close File Handle
close_file <- function(handle) {
  H5Fflush(handle)
  H5Fclose(handle)
}

# Open the Loom file while handling potential locking
open_with_lock <- function(loom_filename, mode) {
  if(mode == "r") mode <- "H5F_ACC_RDONLY"
  if(mode == "r+") mode <- "H5F_ACC_RDWR"
  repeat{ # Handle the lock of the file
    tryCatch({
      return(H5Fopen(name = loom_filename, flags = mode))
    }, error = function(err) {
      if(!grepl("Unable to open file", err$message)) error.json(err$message)
    })
    message("Sleeping 1sec for file lock....")
    time_idle <<- time_idle + 1
    Sys.sleep(1)
  }
}

# Delete a dataset if exists
delete_dataset <- function(handle, dataset_path) {
  tryCatch({
    h5delete(file = handle, name = dataset_path)
  }, error = function(err) {
    if(!grepl("Specified link doesn't exist.", err$message)) error.json(err$message)
  })
}

add_matrix_dataset <- function(handle, dataset_path, dataset_object, storage.mode_param="double", chunk_param=c(min(dim(dataset_object)[1], 64), min(dim(dataset_object)[2], 64)), level_param = 2) {
  if(is.null(dim(dataset_object))) error.json("Cannot write this dataset, it is not a matrix!")
  delete_dataset(handle, dataset_path)
  out <- h5createDataset(file = handle, dataset = dataset_path, dims = dim(dataset_object), storage.mode = storage.mode_param, chunk=chunk_param, level=level_param)
  h5write(file = handle, obj = dataset_object, name=dataset_path)
}

add_array_dataset <- function(handle, dataset_path, dataset_object, storage.mode_param="double", chunk_param=c(min(length(dataset_object), 64)), level_param = 2) {
  if(!is.null(dim(dataset_object))) error.json("Cannot write this dataset, it is not an array!")
  delete_dataset(handle, dataset_path)
  out <- h5createDataset(file = handle, dataset = dataset_path, dims = length(dataset_object), storage.mode = storage.mode_param, chunk=chunk_param, level=level_param)
  h5write(file = handle, obj = dataset_object, name=dataset_path)
}

# Fetch dataset if exists or return NULL
fetch_dataset <- function(handle, dataset_path, transpose = F) {
  handle_dataset <- NULL
  tryCatch({
	handle_dataset <- H5Dopen(h5loc = handle, name = dataset_path)
	out <- H5Dread(handle_dataset)
	H5Dclose(handle_dataset)
    if(transpose) out <- t(out)
    if(length(dim(out)) == 2) return(as.data.frame(out))
    return(as.vector(out))
  }, error = function(err) {
	if(!is.null(handle_dataset)) H5Dclose(handle_dataset)
    if(grepl("Can't open object", err$message)) return(NULL)
    else error.json(err$message)
  })
}
