library(haven)
library(dplyr)
library(cmdstanr)

nSims <- 10000

# running sims in batches
nBatches <- 10

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Grab the batch argument (first argument passed in)
batch <- as.numeric(args[1])

# running sims in nBatches batches
datasetsToAnalyse <- ((batch-1)*(nSims/nBatches)+1):((batch)*(nSims/nBatches))

vars <- c(
  "ageAtOnsetIntOrig", "ageAtOnsetParent", "sigma_ageAtOnset",
  "alpha_0Orig", "alpha_tOrig", "sdYResOrig", "fitted_mean"
)
varsAgeAdj <- c(vars, "alpha_oOrig", "alpha_t_minus_o")

y_crossModelRes <- list(ests = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 stderrs = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 ciLower = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 ciUpper = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 essBulk = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 essTail = array(0, dim=c(length(datasetsToAnalyse),length(vars))))
y_crossModelAgeOnsetAdjRes <- list(ests = array(0, dim=c(length(datasetsToAnalyse),length(varsAgeAdj))),
                                   stderrs = array(0, dim=c(length(datasetsToAnalyse),length(varsAgeAdj))),
                                   ciLower = array(0, dim=c(length(datasetsToAnalyse),length(varsAgeAdj))),
                                   ciUpper = array(0, dim=c(length(datasetsToAnalyse),length(varsAgeAdj))),
                                   essBulk = array(0, dim=c(length(datasetsToAnalyse),length(varsAgeAdj))),
                                   essTail = array(0, dim=c(length(datasetsToAnalyse),length(varsAgeAdj))))

y2_crossModelRes <- y_crossModelRes
y2_crossModelAgeOnsetAdjRes <- y_crossModelAgeOnsetAdjRes

init_fun <- function() {
  list(
    onset_cens = stan_data$censor_age + 1
  )
}

crossModel <- cmdstan_model(exe_file="crossModel.exe")
crossModelAgeOnsetAdj <- cmdstan_model(exe_file="crossModelAgeOnsetAdj.exe")

metrics <- c("ests", "stderrs", "ciLower", "ciUpper", "essBulk", "essTail")


run_crossModel_fit <- function(stan_model, stan_data, init_fun, vars,
                               chains = 1, parallel_chains = 1,
                               iter_warmup = 1000, iter_sampling = 10000,
                               seed = NULL, refresh = 100) {
  
  fit <- stan_model$sample(
    data = stan_data,
    init = init_fun,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh
  )
  
  draws <- fit$draws(variables = vars)
  s <- posterior::summarise_draws(
    draws, mean, sd,
    ~quantile(.x, 0.025), ~quantile(.x, 0.975),
    posterior::ess_bulk, posterior::ess_tail
  )
  
  list(
    variable  = s$variable,
    ests      = s$mean,
    stderrs   = s$sd,
    ciLower   = s$`2.5%`,
    ciUpper   = s$`97.5%`,
    essBulk   = s$`posterior::ess_bulk`,
    essTail   = s$`posterior::ess_tail`
  )
}


j <- 0

for (i in datasetsToAnalyse) {
  j <- j+1
  print(i)
  
  # Stan

  fileName <- paste("Datasets/D", i, ".dta", sep = "")
  simData <- read_stata(fileName)

  stan_data <- list(
    n_subj = nrow(simData),
    
    parentOnset = as.numeric(simData$parental_onset-40),
    
    age = as.numeric(simData$age-40),
    y = (as.numeric(simData$y) - 100) / 10,
    
    n_obs_onset = sum(!is.na(simData$actual_onset_observed)),
    n_cens_onset = sum(is.na(simData$actual_onset_observed)),
    
    subj_obs_onset = as.integer(simData$id[!is.na(simData$actual_onset_observed)]),
    subj_cens_onset = as.integer(simData$id[is.na(simData$actual_onset_observed)]),
    
    onset_obs = as.numeric(simData$actual_onset_observed[!is.na(simData$actual_onset_observed)]-40),
    censor_age = as.numeric(simData$age[is.na(simData$actual_onset_observed)]-40)
  )
  
  # fit models to y
    res <- run_crossModel_fit(crossModel, stan_data, init_fun, vars, seed = i)
    for (m in metrics) y_crossModelRes[[m]][j, ] <- res[[m]]
  
    # fit model with age at onset adjustment
    res <- run_crossModel_fit(crossModelAgeOnsetAdj, stan_data, init_fun, varsAgeAdj, seed = i)
    for (m in metrics) y_crossModelAgeOnsetAdjRes[[m]][j, ] <- res[[m]]
  
  # same models for y2 now
    stan_data$y <- (as.numeric(simData$y2) - 100) / 10
    
    res <- run_crossModel_fit(crossModel, stan_data, init_fun, vars, seed = i)
    for (m in metrics) y2_crossModelRes[[m]][j, ] <- res[[m]]
    
    res <- run_crossModel_fit(crossModelAgeOnsetAdj, stan_data, init_fun, varsAgeAdj, seed = i)
    for (m in metrics) y2_crossModelAgeOnsetAdjRes[[m]][j, ] <- res[[m]]
  
}
  
save(y_crossModelRes,
     y_crossModelAgeOnsetAdjRes,
     y2_crossModelRes,
     y2_crossModelAgeOnsetAdjRes,
       file=paste("cross_",batch,".RData",sep=""))