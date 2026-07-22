library(haven)
library(dplyr)
library(cmdstanr)

nSims <- 1000

# running sims in batches
nBatches <- 500

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Grab the batch argument (first argument passed in)
batch <- as.numeric(args[1])

# running sims in nBatches batches
datasetsToAnalyse <- ((batch-1)*(nSims/nBatches)+1):((batch)*(nSims/nBatches))

vars <- c(
  "ageAtOnsetInt", "sigma_ageAtOnset", "sd_aoo_mutation",
  "alpha0", "b", "sd_alpha[1]", "sd_alpha[2]", "Cor_alpha[1,2]",
  "sd_mutation", "sdYRes"
)

# spline model: parameter names (b is now a length-3 vector; also grab the
# back-transformed original-scale coefficients)
vars_spline <- c(
  "ageAtOnsetInt", "sigma_ageAtOnset", "sd_aoo_mutation",
  "alpha0_orig", "b_orig[1]", "b_orig[2]", "b_orig[3]",
  "sd_alpha[1]", "sd_alpha[2]", "Cor_alpha[1,2]",
  "sd_mutation", "sdYRes"
)

make_res_container <- function(n_reps, n_var) {
  list(ests = array(0, dim=c(n_reps, n_var)),
       stderrs = array(0, dim=c(n_reps, n_var)),
       ciLower = array(0, dim=c(n_reps, n_var)),
       ciUpper = array(0, dim=c(n_reps, n_var)),
       essBulk = array(0, dim=c(n_reps, n_var)),
       essTail = array(0, dim=c(n_reps, n_var)))
}

n_reps <- length(datasetsToAnalyse)

# y1: linear model
y1_Res <- make_res_container(n_reps, length(vars))

# y2: spline model (coefficients + grid predictions)
y2_Spline_Res      <- make_res_container(n_reps, length(vars_spline))
y2_Spline_Pred_Res <- NULL  # sized once n_grid is known, below

# y3: spline model (coefficients + grid predictions)
y3_Spline_Res      <- make_res_container(n_reps, length(vars_spline))
y3_Spline_Pred_Res <- NULL  # sized once n_grid is known, below

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

# spline model: same init logic, but b is a length-3 zero vector, and this
# reads from stan_data_spline (reassigned before each call, for y2 then y3)
init_fun_spline <- function() {
  list(
    ageAtOnsetInt = mean(stan_data_spline$onset_obs),
    sigma_ageAtOnset = max(sd(stan_data_spline$onset_obs), 1),
    
    sd_alpha = c(1, 0.5),
    Lcorr = diag(2),
    z_alpha = matrix(0, nrow = 2, ncol = stan_data_spline$n_subj),
    
    sd_mutation = 1,
    sd_aoo_mutation = 1,
    z_mut = rep(0, stan_data_spline$n_mutation),
    z_aoo_mut = rep(0, stan_data_spline$n_mutation),
    
    alpha0 = mean(stan_data_spline$y),
    b = rep(0, 3),
    sdYRes = max(sd(stan_data_spline$y), 0.5),
    
    onset_cens = censor_age + 1
  )
}

longModelLinear <- cmdstan_model(exe_file="longModelLinear")

longModelSpline <- cmdstan_model(exe_file="longModelSpline")

# preparations for cubic spline model
knot <- c(-20,-10,0,10)

# pre-calculate QR decomposition based on parental age of onset in first simulated dataset
fileName <- paste("Datasets/NewData1.dta",sep="")
simData <- read_stata(fileName)
# only keep rows where age<= final_visit_age and generation==1
simData <- subset(simData, age<= final_visit_age)
simData <- subset(simData, generation==1)

timeSinceParental <- simData$age-simData$parental_onset

term1Pre <- pmax((timeSinceParental-knot[1])^3,0)
term2Pre <- pmax((timeSinceParental-knot[2])^3,0)
term3Pre <- pmax((timeSinceParental-knot[3])^3,0)
term4Pre <- pmax((timeSinceParental-knot[4])^3,0)

u1Pre <- timeSinceParental
u2Pre = (term1Pre-(term3Pre*(knot[4]-knot[1])-term4Pre*(knot[3]-knot[1]))/(knot[4]-knot[3]))/(knot[4]-knot[1])^2
u3Pre = (term2Pre-(term3Pre*(knot[4]-knot[2])-term4Pre*(knot[3]-knot[2]))/(knot[4]-knot[3]))/(knot[4]-knot[1])^2

uMat <- cbind(u1Pre, u2Pre, u3Pre)
cor(uMat)
# standardize to have mean zero and SD 1
u_means <- colMeans(uMat)
u_sds <- apply(uMat, 2, sd)
u_tilde <- scale(uMat, center = u_means, scale = u_sds)

qr_u <- qr(u_tilde)
Q <- qr.Q(qr_u)  # n x p matrix (thin Q)
cor(Q)
R <- qr.R(qr_u)  # p x p upper-triangular matrix

# Transformation matrix: T = inv(R).
T <- solve(R)
#Ws <- u_tilde %*% T
#print(all.equal(Ws, Q)) 

# set up matrix for getting predictions at grid of time relative to onset
timeGrid <- seq(-30,30,1)
term1Grid <- pmax((timeGrid-knot[1])^3,0)
term2Grid <- pmax((timeGrid-knot[2])^3,0)
term3Grid <- pmax((timeGrid-knot[3])^3,0)
term4Grid <- pmax((timeGrid-knot[4])^3,0)
u1Grid <- timeGrid
u2Grid = (term1Grid-(term3Grid*(knot[4]-knot[1])-term4Grid*(knot[3]-knot[1]))/(knot[4]-knot[3]))/(knot[4]-knot[1])^2
u3Grid = (term2Grid-(term3Grid*(knot[4]-knot[2])-term4Grid*(knot[3]-knot[2]))/(knot[4]-knot[3]))/(knot[4]-knot[1])^2
uMatGrid <- cbind(u1Grid, u2Grid, u3Grid)
u_tildeGrid <- scale(uMatGrid, center = u_means, scale = u_sds)
gridDesign <- u_tildeGrid %*% T

# strip attributes and rename to match the Stan data block's expected field
# names ("splineRot" for T, "t_grid"/"n_grid" for the grid)
splineRot <- matrix(as.numeric(T), nrow = 3, ncol = 3)
t_grid <- as.numeric(timeGrid)
n_grid <- length(t_grid)

vars_pred <- paste0("y_pred[", 1:n_grid, "]")

# now size the grid prediction containers, since n_grid is known
y2_Spline_Pred_Res <- make_res_container(n_reps, n_grid)
y3_Spline_Pred_Res <- make_res_container(n_reps, n_grid)

metrics <- c("ests", "stderrs", "ciLower", "ciUpper", "essBulk", "essTail")

run_longModel_fit <- function(stan_model, stan_data, init_fun, vars,
                              chains = 1, parallel_chains = 1,
                              iter_warmup = 1000, iter_sampling = 1000,
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
    variable = s$variable,
    ests = s$mean,
    stderrs = s$sd,
    ciLower = s$`2.5%`,
    ciUpper = s$`97.5%`,
    essBulk = s$`posterior::ess_bulk`,
    essTail = s$`posterior::ess_tail`,
    fit = fit   # keep the fit object so callers can pull additional variables without refitting
  )
}

# fits the spline model for a given outcome column, storing coefficient and
# grid-prediction results into the two supplied containers at row j
fit_spline_outcome <- function(base_stan_data, y_col, coef_res, pred_res, j, seed) {
  stan_data_spline <<- c(base_stan_data, list(
    y = (as.numeric(y_col) - 100) / 10,
    knot = as.numeric(knot),
    u_means = as.numeric(u_means),
    u_sds = as.numeric(u_sds),
    splineRot = splineRot,
    t_grid = t_grid,
    n_grid = n_grid
  ))
  
  res_spline <- run_longModel_fit(longModelSpline, stan_data_spline, init_fun_spline, vars_spline, seed=seed)
  for (m in metrics) coef_res[[m]][j, ] <- res_spline[[m]]
  
  draws_pred <- res_spline$fit$draws(variables = vars_pred)
  s_pred <- posterior::summarise_draws(
    draws_pred, mean, sd,
    ~quantile(.x, 0.025), ~quantile(.x, 0.975),
    posterior::ess_bulk, posterior::ess_tail
  )
  pred_res$ests[j, ]    <- s_pred$mean
  pred_res$stderrs[j, ] <- s_pred$sd
  pred_res$ciLower[j, ] <- s_pred$`2.5%`
  pred_res$ciUpper[j, ] <- s_pred$`97.5%`
  pred_res$essBulk[j, ] <- s_pred$`posterior::ess_bulk`
  pred_res$essTail[j, ] <- s_pred$`posterior::ess_tail`
  
  list(coef_res = coef_res, pred_res = pred_res, fit = res_spline$fit)
}

j <- 0

for (i in datasetsToAnalyse) {
  j <- j+1
  print(i)
  
  # load sim data
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
  
  # fields common to all three outcome models (everything except y itself)
  base_stan_data <- list(
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
    
    n_obs_onset = length(subj_obs_onset),
    n_cens_onset = length(subj_cens_onset),
    
    subj_obs_onset = as.integer(subj_obs_onset),
    subj_cens_onset = as.integer(subj_cens_onset),
    
    onset_obs = as.numeric(onset_obs),
    censor_age = as.numeric(censor_age)
  )
  
  # --- y1: linear model ---
  stan_data <- c(base_stan_data, list(y = (as.numeric(simData$y1) - 100) / 10))
  res <- run_longModel_fit(longModelLinear, stan_data, init_fun, vars, seed=i)
  for (m in metrics) y1_Res[[m]][j, ] <- res[[m]]
  
  # --- y2: spline model ---
  out2 <- fit_spline_outcome(base_stan_data, simData$y2, y2_Spline_Res, y2_Spline_Pred_Res, j, seed=i)
  y2_Spline_Res      <- out2$coef_res
  y2_Spline_Pred_Res <- out2$pred_res
  
  # --- y3: spline model ---
  out3 <- fit_spline_outcome(base_stan_data, simData$y3, y3_Spline_Res, y3_Spline_Pred_Res, j, seed=i)
  y3_Spline_Res      <- out3$coef_res
  y3_Spline_Pred_Res <- out3$pred_res
  
}

save(y1_Res,
     y2_Spline_Res, y2_Spline_Pred_Res,
     y3_Spline_Res, y3_Spline_Pred_Res,
     file=paste("long_",batch,".RData",sep=""))
