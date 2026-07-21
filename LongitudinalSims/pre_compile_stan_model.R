library(cmdstanr)

longModelLinear <- cmdstan_model("longModelLinear.stan")
longModelSpline <- cmdstan_model("longModelSpline.stan")
