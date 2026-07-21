
nBatches <- 10
prefix <- "long"

batchLoad <- new.env()

# load first batch into global env
load(paste0(prefix, "_1.RData"))

listNames <- c("y1_Res")
  
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
  46,
  9*sqrt(0.4),
  9*sqrt(0.6),
  0,
  0.1,
  sqrt(5^2+25^2*0.5^2)/10,
  0.5/10,
  (25*0.5^2)/(sqrt(5^2+5^2+25^2*0.5^2)*sqrt(0.5^2)),
  0.5,
  0.5)

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

summarise_res(y1_Res,  trueParms_set1)


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
