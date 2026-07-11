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
  "alpha0Orig", "bOrig", "sdYResOrig", "fitted_mean"
)
vars2 <- c(
  "ageAtOnsetIntOrig", "ageAtOnsetParent", "sigma_ageAtOnset",
  "alpha0Orig", "bOrig", "cOrig", "sdYResOrig", "c_minus_b", "fitted_mean"
)

yMod1Res <- list(ests = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 stderrs = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 ciLower = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 ciUpper = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 essBulk = array(0, dim=c(length(datasetsToAnalyse),length(vars))),
                 essTail = array(0, dim=c(length(datasetsToAnalyse),length(vars))))
y2Mod1Res <- yMod1Res
yMod2Res <- list(ests = array(0, dim=c(length(datasetsToAnalyse),length(vars2))),
                 stderrs = array(0, dim=c(length(datasetsToAnalyse),length(vars2))),
                 ciLower = array(0, dim=c(length(datasetsToAnalyse),length(vars2))),
                 ciUpper = array(0, dim=c(length(datasetsToAnalyse),length(vars2))),
                 essBulk = array(0, dim=c(length(datasetsToAnalyse),length(vars2))),
                 essTail = array(0, dim=c(length(datasetsToAnalyse),length(vars2))))
y2Mod2Res <- yMod2Res


init_fun <- function() {
  list(
    ageAtOnsetInt = 0,
    ageAtOnsetParent = 0,
    sigma_ageAtOnset = max(sd(stan_data$onset_obs), 1),
    
    alpha0 = mean(stan_data$y),
    b = 0,
    sdYRes = max(sd(stan_data$y), 0.5),
    
    onset_cens = stan_data$censor_age + 1
  )
}
init_fun2 <- function() {
  list(
    ageAtOnsetInt = 0,
    ageAtOnsetParent = 0,
    sigma_ageAtOnset = max(sd(stan_data$onset_obs), 1),
    
    alpha0 = mean(stan_data$y),
    b = 0,
    c = 0,
    sdYRes = max(sd(stan_data$y), 0.5),
    
    onset_cens = stan_data$censor_age + 1
  )
}

mod1 <- cmdstan_model("crossModel1.stan")
mod2 <- cmdstan_model("crossModel2.stan")

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
  
  # model 1 to y
  
  fit <- mod1$sample(
    data = stan_data,
    init = init_fun,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 1000,
    iter_sampling = 10000,
    seed = i,
    refresh = 100
  )
  
  # fit$diagnostic_summary()

  draws <- fit$draws(variables = vars)
  s <- posterior::summarise_draws(draws, mean, sd, ~quantile(.x, 0.025), ~quantile(.x, 0.975),
                                  posterior::ess_bulk,
                                  posterior::ess_tail)
  
  yMod1Res$ests[j, ] <- s$mean
  yMod1Res$stderrs[j, ] <- s$sd
  yMod1Res$ciLower[j, ] <- s$`2.5%`
  yMod1Res$ciUpper[j, ] <- s$`97.5%`
  yMod1Res$essBulk[j, ] <- s$`posterior::ess_bulk`
  yMod1Res$essTail[j, ] <- s$`posterior::ess_tail`

  # y model 2
  
  fit <- mod2$sample(
    data = stan_data,
    init = init_fun2,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 1000,
    iter_sampling = 10000,
    seed = i,
    refresh = 100
  )

  draws <- fit$draws(variables = vars2)
  s <- posterior::summarise_draws(draws, mean, sd, ~quantile(.x, 0.025), ~quantile(.x, 0.975),
                                  posterior::ess_bulk,
                                  posterior::ess_tail)

  yMod2Res$ests[j, ] <- s$mean
  yMod2Res$stderrs[j, ] <- s$sd
  yMod2Res$ciLower[j, ] <- s$`2.5%`
  yMod2Res$ciUpper[j, ] <- s$`97.5%`
  yMod2Res$essBulk[j, ] <- s$`posterior::ess_bulk`
  yMod2Res$essTail[j, ] <- s$`posterior::ess_tail`
  
  
  # y2 model 1
  
  stan_data$y <- (as.numeric(simData$y2) - 100) / 10
  
  fit <- mod1$sample(
    data = stan_data,
    init = init_fun,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 1000,
    iter_sampling = 10000,
    seed = i,
    refresh = 100
  )
  
  draws <- fit$draws(variables = vars)
  s <- posterior::summarise_draws(draws, mean, sd, ~quantile(.x, 0.025), ~quantile(.x, 0.975),
                                  posterior::ess_bulk,
                                  posterior::ess_tail)
  
  y2Mod1Res$ests[j, ] <- s$mean
  y2Mod1Res$stderrs[j, ] <- s$sd
  y2Mod1Res$ciLower[j, ] <- s$`2.5%`
  y2Mod1Res$ciUpper[j, ] <- s$`97.5%`
  y2Mod1Res$essBulk[j, ] <- s$`posterior::ess_bulk`
  y2Mod1Res$essTail[j, ] <- s$`posterior::ess_tail`
  
  # y2 model 2
  
  fit <- mod2$sample(
    data = stan_data,
    init = init_fun2,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 1000,
    iter_sampling = 10000,
    seed = i,
    refresh = 100
  )
  
  draws <- fit$draws(variables = vars2)
  s <- posterior::summarise_draws(draws, mean, sd, ~quantile(.x, 0.025), ~quantile(.x, 0.975),
                                  posterior::ess_bulk,
                                  posterior::ess_tail)
  
  y2Mod2Res$ests[j, ] <- s$mean
  y2Mod2Res$stderrs[j, ] <- s$sd
  y2Mod2Res$ciLower[j, ] <- s$`2.5%`
  y2Mod2Res$ciUpper[j, ] <- s$`97.5%`
  y2Mod2Res$essBulk[j, ] <- s$`posterior::ess_bulk`
  y2Mod2Res$essTail[j, ] <- s$`posterior::ess_tail`
  
}
  
save(yMod1Res,yMod2Res,y2Mod1Res,y2Mod2Res,
       file=paste("cross_",batch,".RData",sep=""))