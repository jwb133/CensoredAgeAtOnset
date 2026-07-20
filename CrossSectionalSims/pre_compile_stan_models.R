library(cmdstanr)

crossModel <- cmdstan_model("crossModel.stan")
crossModelAgeAdj <- cmdstan_model("crossModelAgeAdj.stan")
crossModelQuad <- cmdstan_model("crossModelQuad.stan")
