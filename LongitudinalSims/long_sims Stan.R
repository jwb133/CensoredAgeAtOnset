library(haven)
library(dplyr)
library(cmdstanr)

nSims <- 1000

# running sims in batches
nBatches <- 10

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Grab the batch argument (first argument passed in)
batch <- as.numeric(args[1])



# running sims in nBatches batches
datasetsToAnalyse <- ((batch-1)*(nSims/nBatches)+1):((batch)*(nSims/nBatches))

ests_linear <- array(0, dim=c(length(datasetsToAnalyse),10))
stderrs_linear <- ests_linear
ciLower_linear <- ests_linear
ciUpper_linear <- ests_linear
essBulk_linear <- ests_linear
essTail_linear <- ests_linear

vars <- c(
  "ageAtOnsetInt", "sigma_ageAtOnset", "sd_aoo_mutation",
  "alpha0", "b", "sd_alpha[1]", "sd_alpha[2]", "Cor_alpha[1,2]",
  "sd_mutation", "sdYRes"
)


init_fun <- function() {
  list(
    ageAtOnsetInt = mean(stan_data$onset_obs),
    sigma_ageAtOnset = max(sd(stan_data$onset_obs), 1),
    
    sd_alpha = c(1, 0.5),
    Lcorr = diag(2),
    z_alpha = matrix(0, nrow = 2, ncol = stan_data$n_subj),
    
    sd_mutation = 1,
    sd_aoo_mutation = 1,
    z_mut = rep(0, stan_data$n_mutation),
    z_aoo_mut = rep(0, stan_data$n_mutation),
    
    alpha0 = mean(stan_data$y),
    b = 0,
    sdYRes = max(sd(stan_data$y), 0.5),
    
    onset_cens = censor_age + 1
  )
}

mod <- cmdstan_model("longModelLinear.stan")

j <- 0

for (i in datasetsToAnalyse) {
  j <- j+1
  print(i)
  
  # Stan

  fileName <- paste("Datasets/NewData", i, ".dta", sep = "")
  simData <- read_stata(fileName)
  
  # create new mutation ids
  simData <- simData %>% mutate(mutation = dense_rank(mutation))
  
  # parent onset data: generation 0
  parentOnsetData <- simData[simData$generation == 0, ]
  n_mutation <- max(simData$mutation)
  
  # keep visits up to final visit age
  simData <- subset(simData, age <= final_visit_age)
  
  # keep only generation 1 carriers for longitudinal model
  simData <- subset(simData, generation == 1)
  simData <- subset(simData, carrier == 1)
  
  # create new subject ids starting at 1
  simData <- simData %>% mutate(id = dense_rank(id))
  
  # subject-level dataset
  timeInvariant <- data.frame(
    id = simData$id,
    actual_onset_observed = simData$actual_onset_observed,
    final_visit_age = simData$final_visit_age,
    mutation = simData$mutation
  )
  
  timeInvariant <- unique(timeInvariant)
  timeInvariant <- timeInvariant[order(timeInvariant$id), ]
  
  # observed/censored split
  onset_obs_full <- as.numeric(timeInvariant$actual_onset_observed)
  final_visit_age_full <- as.numeric(timeInvariant$final_visit_age)
  
  # observed if actual_onset_observed is not NA
  is_obs <- !is.na(onset_obs_full)
  
  subj_obs_onset  <- which(is_obs)
  subj_cens_onset <- which(!is_obs)
  
  onset_obs  <- onset_obs_full[subj_obs_onset]
  censor_age <- final_visit_age_full[subj_cens_onset]
  
  # checks
  stopifnot(length(subj_obs_onset) + length(subj_cens_onset) == nrow(timeInvariant))
  stopifnot(all(!is.na(onset_obs)))
  stopifnot(all(!is.na(censor_age)))
  stopifnot(all(sort(unique(simData$id)) == seq_len(nrow(timeInvariant))))
  
  # optional stronger check:
  # censored subjects should have no observed onset
  stopifnot(all(is.na(onset_obs_full[subj_cens_onset])))
  
  stan_data <- list(
    n_subj = nrow(timeInvariant),
    n_obs = nrow(simData),
    n_mutation = n_mutation,
    n_parents = nrow(parentOnsetData),
    
    parentOnset = as.numeric(parentOnsetData$actual_onset),
    
    mutation_by_sub = as.integer(timeInvariant$mutation),
    mutation_by_parent = as.integer(parentOnsetData$mutation),
    id = as.integer(simData$id),
    mutation = as.integer(simData$mutation),
    
    age = as.numeric(simData$age),
    y = (as.numeric(simData$y1) - 100) / 10,
    
    n_obs_onset = length(subj_obs_onset),
    n_cens_onset = length(subj_cens_onset),
    
    subj_obs_onset = as.integer(subj_obs_onset),
    subj_cens_onset = as.integer(subj_cens_onset),
    
    onset_obs = as.numeric(onset_obs),
    censor_age = as.numeric(censor_age)
  )
  
  fit <- mod$sample(
    data = stan_data,
    init = init_fun,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    seed = i,
    refresh = 100
  )
  
  fit$diagnostic_summary()
  
  
  s <- fit$summary(variables = vars)
  s[, c("variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")]
  
  draws <- fit$draws(variables = vars)
  s <- posterior::summarise_draws(draws, mean, sd, ~quantile(.x, 0.025), ~quantile(.x, 0.975),
                                  posterior::ess_bulk,
                                  posterior::ess_tail)
  
  ests_linear[j, ]    <- s$mean
  stderrs_linear[j, ] <- s$sd
  ciLower_linear[j, ] <- s$`2.5%`
  ciUpper_linear[j, ] <- s$`97.5%`
  essBulk_linear[j, ] <- s$`posterior::ess_bulk`
  essTail_linear[j, ] <- s$`posterior::ess_tail`
  
}
  
save(ests_linear,stderrs_linear,ciLower_linear,ciUpper_linear,
       file=paste("long_",batch,".RData",sep=""))