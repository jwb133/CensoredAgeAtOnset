
nBatches <- 1000
prefix <- "long"

batchLoad <- new.env()

# load first batch into global env
load(paste0(prefix, "_1.RData"))

listNames <- c("y1_Res","y2_Spline_Res","y2_Spline_Pred_Res",
	"y3_Spline_Res","y3_Spline_Pred_Res")
  
components <- c("ests", "stderrs", "ciLower", "ciUpper", "essBulk", "essTail")

for (i in 2:nBatches) {
  load(paste0(prefix, "_", i, ".RData"), envir = batchLoad)
  
  for (ln in listNames) {
    for (comp in components) {
      
      # current (accumulated) data
      current <- get(ln, envir = .GlobalEnv)[[comp]]
      
      # new batch data
      new <- get(ln, envir = batchLoad)[[comp]]
      
      # rbind and store back
      tmp <- get(ln, envir = .GlobalEnv)
      tmp[[comp]] <- rbind(current, new)
      assign(ln, tmp, envir = .GlobalEnv)
    }
  }
}

save(y1_Res,
     y2_Spline_Res, y2_Spline_Pred_Res,
     y3_Spline_Res, y3_Spline_Pred_Res,
     file = paste0(prefix, "_combined.RData"))

for (obj_name in listNames) {
  # Get the object from the global environment
  mat <- get(obj_name, envir = .GlobalEnv)
  mat$essBulk <- NULL
  mat$essTail <- NULL
  
  # Construct filename
  filename <- paste0(prefix, "_", obj_name, ".csv")
  
  # Write to CSV
  write.csv(mat, file = filename, row.names = FALSE)
}
