
nBatches <- 10
prefix <- "cross"

batchLoad <- new.env()

# load first batch into global env
load(paste0(prefix, "_1.RData"))

listNames <- c("y_crossModelRes",
               "y_crossModelAgeOnsetAdjRes",
               "y2_crossModelRes",
               "y2_crossModelAgeOnsetAdjRes")
  
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

# evaluate results

# true parameters in model 1 when model 1 is correct (no age at onset effect)
trueParms_set1 <- c(
  46 * (1 - 48.6 / 81),
  48.6 / 81,
  sqrt(81 * (1 - (48.6 / 81)^2)),
  100,
  1,
  27,
  94
)

trueParms_set2 <- c(
  46 * (1 - 48.6 / 81),
  48.6 / 81,
  sqrt(81 * (1 - (48.6 / 81)^2)),
  100,
  1,
  27,
  94,
  0,
  1
)

trueParms_set4 <- c(
  46 * (1 - 48.6 / 81),
  48.6 / 81,
  sqrt(81 * (1 - (48.6 / 81)^2)),
  100,
  2,
  27,
  94,
  1,
  1
)



summarise_res <- function(res, trueParms) {
  list(
    mean_est = colMeans(res$ests),
    true = trueParms,
    emp_sd = apply(res$ests, 2, sd),
    mean_post_sd = colMeans(res$stderrs),
    coverage = sapply(seq_along(trueParms), function(i) {
      mean((res$ciLower[, i] < trueParms[i]) & (res$ciUpper[, i] > trueParms[i]))
    }),
    mean_ess_bulk = colMeans(res$essBulk),
    mean_ess_tail = colMeans(res$essTail)
  )
}

summarise_res(y_crossModelRes,  trueParms_set1)
summarise_res(y_crossModelAgeOnsetAdjRes, trueParms_set2)
summarise_res(y2_crossModelAgeOnsetAdjRes, trueParms_set4)

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
