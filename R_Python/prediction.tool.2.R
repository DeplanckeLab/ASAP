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
  #write(jsonlite::toJSON(stats, method="C", auto_unbox=T, digits = NA), file = paste0(output_dir, "output.json"), append=F)
  message(toJSON(stats, method="C", auto_unbox=T))
  close_all()
  stop()
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
  #message("Reading JSON...")
  if(!file.exists(input_json)) error_json("The input JSON file does not exist!")
  data.json <- read_json(input_json, simplifyVector = TRUE)
  ##message(nrow(data.json), " methods were found.")
  for(i in 1:nrow(data.json)){ # For every method
    data.runs <- data.json[i, "runs"][[1]]
    #message(nrow(data.runs), " runs were found for method '", data.json[i, "std_method_name"], "'")
    
    # 1. Predicting t
    data.runs.t <- data.runs[!is.na(data.runs$t) & !is.na(data.runs$nber_rows) & !is.na(data.runs$nber_cols),]
    #diff <- nrow(data.runs) - nrow(data.runs.t)
    #if(diff > 0) message("[Warning!] ", diff, " runs have NAs and cannot be used for predicting t")
    if(nrow(data.runs.t) < 3) {
      #message("There are less than 3 runs for this method (",nrow(data.runs.t),"), no prediction model can be created")
    } else {
      #message("Building model for '", data.json[i, "std_method_name"], "' (id=", data.json[i, "std_method_id"],") using ", nrow(data.runs.t), " complete runs (without NAs)...")
      lm.model <- lm(data = data.runs.t, formula = t ~ nber_cols * nber_rows)
      model.path <- paste0(model_dir, "model_",data.json[i, "std_method_id"],".t.rda")
      #message("Saving model to ", model.path)
      save(lm.model, file = model.path)
      #predicted_df <- data.frame(t_pred = predict(object = lm.model, data = data.runs.t), dim_both=data.runs.t$nber_cols * data.runs.t$nber_rows, dim_col=data.runs.t$nber_cols, dim_row=data.runs.t$nber_rows)
      #predicted_df = predicted_df[with(predicted_df, order(dim_both)),]
      #ggplot(data.runs.t, aes(x = nber_cols * nber_rows, y = t)) + geom_point() + geom_line(data = predicted_df, mapping = aes(x=dim_both, y=t_pred), color='red')
      #predicted_df = predicted_df[with(predicted_df, order(dim_col)),]
      #ggplot(data.runs.t, aes(x = nber_cols, y = t)) + geom_point() + geom_line(data = predicted_df, mapping = aes(x=dim_col, y=t_pred), color='red')
      #predicted_df = predicted_df[with(predicted_df, order(dim_row)),]
      #ggplot(data.runs.t, aes(x = nber_rows, y = t)) + geom_point() + geom_line(data = predicted_df, mapping = aes(x=dim_row, y=t_pred), color='red')
    }
    
    # 2. Predicting RAM
    data.runs.m <- data.runs[!is.na(data.runs$m) & !is.na(data.runs$nber_rows) & !is.na(data.runs$nber_cols),]
    #diff <- nrow(data.runs) - nrow(data.runs.m)
    #if(diff > 0) message("[Warning!] ", diff, " runs have NAs and cannot be used for predicting RAM")
    if(nrow(data.runs.m) < 3) {
      #message("There are less than 3 runs for this method (",nrow(data.runs.m),"), no prediction model can be created")
    } else {
      #message("Building model for '", data.json[i, "std_method_name"], "' (id=", data.json[i, "std_method_id"],") using ", nrow(data.runs.m), " complete runs (without NAs)...")
      lm.model <- lm(data = data.runs.m, formula = m ~ nber_cols * nber_rows)
      model.path <- paste0(model_dir, "model_",data.json[i, "std_method_id"],".m.rda")
      #message("Saving model to ", model.path)
      save(lm.model, file = model.path)
      #predicted_df <- data.frame(t_pred = predict(object = lm.model, data = data.runs.m), dim_both=data.runs.m$nber_cols * data.runs.m$nber_rows, dim_col=data.runs.m$nber_cols, dim_row=data.runs.m$nber_rows)
      #predicted_df = predicted_df[with(predicted_df, order(dim_both)),]
      #ggplot(data.runs.m, aes(x = nber_cols * nber_rows, y = t)) + geom_point() + geom_line(data = predicted_df, mapping = aes(x=dim_both, y=t_pred), color='red')
      #predicted_df = predicted_df[with(predicted_df, order(dim_col)),]
      #ggplot(data.runs.m, aes(x = nber_cols, y = t)) + geom_point() + geom_line(data = predicted_df, mapping = aes(x=dim_col, y=t_pred), color='red')
      #predicted_df = predicted_df[with(predicted_df, order(dim_row)),]
      #ggplot(data.runs.m, aes(x = nber_rows, y = t)) + geom_point() + geom_line(data = predicted_df, mapping = aes(x=dim_row, y=t_pred), color='red')
    }
  }
} else if (mode == "predict"){
  std_method_id <- as.integer(args[3])
  nber_rows <- as.integer(args[4])
  nber_cols <- as.integer(args[5])
  
  predict_df <- data.frame(nber_rows = nber_rows, nber_cols = nber_cols, m = NA, t = NA)
  
  model.path <- paste0(model_dir, "model_",std_method_id,".t.rda")
  #message("Loading t model from ", model.path)
  predicted_t <- NA
  if(file.exists(model.path)){
    load(file = model.path)
    predicted_t <- NA
    tryCatch ({ predicted_t <<- as.numeric(predict(object = lm.model, newdata = predict_df)) }, warning = function(cond){})
    #message("Predicted time t = ", predicted_t, " s")
  }

  model.path <- paste0(model_dir, "model_",std_method_id,".m.rda")
  #message("Loading m model from ", model.path)
  predicted_m <- NA
  if(file.exists(model.path)){
    load(file = model.path)
    predicted_m <- NA
    tryCatch ({ predicted_m <<- as.numeric(predict(object = lm.model, newdata = predict_df)) }, warning = function(cond){})
    #message("Predicted max RAM m = ", predicted_m, " octets")
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

