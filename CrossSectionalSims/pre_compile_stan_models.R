library(cmdstanr)

crossModel <- cmdstan_model("crossModel.stan")
crossModelAgeAdj <- cmdstan_model("crossModelAgeOnsetAdj.stan")
