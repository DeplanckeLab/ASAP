### ASAP Clustering script
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

### Libraries
suppressPackageStartupMessages(require(jsonlite))
suppressPackageStartupMessages(require(ggplot2))

### Set message() to output in standard output (instead of error output)
sink(stdout(), type = "message")

### Functions
error_json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  message(toJSON(stats, method="C", auto_unbox=T))
  quit(save = "no")
}

### Default Parameters
set.seed(42)
mode <- args[1]
model_dir <- args[2]
if(!endsWith(model_dir, suffix = "/")) model_dir <- paste0(model_dir, "/")
if(!file.exists(model_dir)) error_json("The model folder does not exist!")

### Running mode 1 or 2
if (mode == "build"){
  input_json <- args[3]
  if(!file.exists(input_json)) error_json("The input JSON file does not exist!")
  data.json <- read_json(input_json, simplifyVector = TRUE)
  for(i in 1:nrow(data.json)){ # For every method
    data.runs <- data.json[i, "runs"][[1]]
    data.runs$nber_cols <- as.numeric(data.runs$nber_cols)
    data.runs$nber_rows <- as.numeric(data.runs$nber_rows)
    
    # 1. Predicting t
    data.runs.t <- data.runs[!is.na(data.runs$t) & !is.na(data.runs$nber_rows) & !is.na(data.runs$nber_cols),]
    if(nrow(data.runs.t) >= 3) {
      lm.model <- lm(data = data.runs.t, formula = t ~ nber_cols * nber_rows)
      model.path <- paste0(model_dir, "model_",data.json[i, "std_method_id"],".t.rds")
      saveRDS(lm.model, file = model.path)
    }
    
    # 2. Predicting RAM
    data.runs.m <- data.runs[!is.na(data.runs$m) & !is.na(data.runs$nber_rows) & !is.na(data.runs$nber_cols),]
    if(nrow(data.runs.m) >= 3) {
      lm.model <- lm(data = data.runs.m, formula = m ~ nber_cols * nber_rows)
      model.path <- paste0(model_dir, "model_",data.json[i, "std_method_id"],".m.rds")
      saveRDS(lm.model, file = model.path)
    }
  }
} else if (mode == "predict"){
  std_method_id <- as.integer(args[3])
  nber_rows <- as.integer(args[4])
  nber_cols <- as.integer(args[5])
  
  predict_df <- data.frame(nber_rows = nber_rows, nber_cols = nber_cols, m = NA, t = NA)
  
  model.path <- paste0(model_dir, "model_",std_method_id,".t.rds")
  predicted_t <- NA
  if(file.exists(model.path)){
    lm.model <- readRDS(file = model.path)
    predicted_t <- NA
    tryCatch ({ predicted_t <<- as.numeric(predict(object = lm.model, newdata = predict_df)) }, warning = function(cond){})
  }

  model.path <- paste0(model_dir, "model_",std_method_id,".m.rds")
  predicted_m <- NA
  if(file.exists(model.path)){
    lm.model <- readRDS(file = model.path)
    predicted_m <- NA
    tryCatch ({ predicted_m <<- as.numeric(predict(object = lm.model, newdata = predict_df)) }, warning = function(cond){})
  }
  
  # Generate default JSON file
  if(!is.na(predicted_m) & predicted_m < 0) predicted_m <- 0
  if(!is.na(predicted_t) & predicted_t < 0) predicted_t <- 0
  
  stats <- list()
  stats$predicted_ram <- round(predicted_m)
  stats$predicted_time <- round(predicted_t)
  #write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"), append=F)
  message(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA))
}

