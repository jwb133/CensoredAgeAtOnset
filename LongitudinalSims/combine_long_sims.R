
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

vars <- c(
  "ageAtOnsetInt", "sigma_ageAtOnset", "sd_aoo_mutation",
  "alpha0", "b", "sd_alpha[1]", "sd_alpha[2]", "Cor_alpha[1,2]",
  "sd_mutation", "sdYRes"
)

vars_spline <- c(
  "ageAtOnsetInt", "sigma_ageAtOnset", "sd_aoo_mutation",
  "alpha0_orig", "b_orig[1]", "b_orig[2]", "b_orig[3]",
  "sd_alpha[1]", "sd_alpha[2]", "Cor_alpha[1,2]",
  "sd_mutation", "sdYRes"
)

varSetMap <- list(
  y1_Res             = vars,
  y2_Spline_Res      = vars_spline,
  y2_Spline_Pred_Res = vars_spline,
  y3_Spline_Res      = vars_spline,
  y3_Spline_Pred_Res = vars_spline
)

predObjNames <- c("y2_Spline_Pred_Res", "y3_Spline_Pred_Res")

for (obj_name in listNames) {
  # Get the object from the global environment
  mat <- get(obj_name, envir = .GlobalEnv)
  mat$essBulk <- NULL
  mat$essTail <- NULL

  if (obj_name %in% predObjNames) {
    # Prediction grids: just number the columns per component
    df <- do.call(cbind, lapply(names(mat), function(comp) {
      m <- mat[[comp]]
      colnames(m) <- paste0(comp, "_", seq_len(ncol(m)))
      m
    }))
  } else {
    # Coefficient matrices: use named variables
    varNames <- varSetMap[[obj_name]]
    df <- do.call(cbind, lapply(names(mat), function(comp) {
      m <- mat[[comp]]
      colnames(m) <- paste0(comp, "_", varNames)
      m
    }))
  }

  # Construct filename
  filename <- paste0(prefix, "_", obj_name, ".csv")

  # Write to CSV
  write.csv(df, file = filename, row.names = FALSE)
}
