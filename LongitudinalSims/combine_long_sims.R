prefix <- "long"

# running sims in batches
nBatches <- 10

batchLoad <- new.env()

load(paste0(prefix, "_1.RData"))

wanted <- c("ests_linear", "stderrs_linear", "ciLower_linear", "ciUpper_linear",
            "essBulk_linear","essTail_linear")

for (i in 2:nBatches) {
  load(paste0(prefix, "_",i, ".RData"), envir = batchLoad)
  
  for (obj in intersect(ls(batchLoad), wanted)) {
    assign(
      obj,
      rbind(get(obj, envir = .GlobalEnv),
            get(obj, envir = batchLoad)),
      envir = .GlobalEnv
    )
  }
}

for (obj_name in wanted) {
  # Get the object from the global environment
  mat <- get(obj_name, envir = .GlobalEnv)
  
  # Construct filename
  filename <- paste0(prefix, "_", obj_name, ".csv")
  
  # Write to CSV
  write.csv(mat, file = filename, row.names = FALSE)
}

write.csv(linear_parms, file=paste0(prefix,"_linear_parms.csv"))
write.csv(spline_parms, file=paste0(prefix,"_spline_parms.csv"))