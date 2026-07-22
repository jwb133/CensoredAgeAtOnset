# load results datasets
load("long_combined.RData")

# evaluate simulations

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

trueParms_set2 <- c(
  46,
  9*sqrt(0.4),
  9*sqrt(0.6),
  0,
  0,
  0,
  0.15,
  sqrt(5^2+25^2*0.5^2)/10,
  0.5/10,
  (25*0.5^2)/(sqrt(5^2+5^2+25^2*0.5^2)*sqrt(0.5^2)),
  0.5,
  0.5)

trueParms_set3 <- c(
  46,
  9*sqrt(0.4),
  9*sqrt(0.6),
  0,
  0.1,
  0,
  0.15,
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

summarise_res(y1_Res, trueParms_set1)
summarise_res(y2_Spline_Res, trueParms_set2)
summarise_res(y3_Spline_Res, trueParms_set3)

# evaluate spline fits at grid
# first re-generate true values for the three y's across grid
knot <- c(-20, -10, 0, 10)
timeGrid <- seq(-30, 30, 1)

term1Grid <- pmax((timeGrid - knot[1])^3, 0)
term2Grid <- pmax((timeGrid - knot[2])^3, 0)
term3Grid <- pmax((timeGrid - knot[3])^3, 0)
term4Grid <- pmax((timeGrid - knot[4])^3, 0)

u1Grid <- timeGrid
u2Grid <- (term1Grid - (term3Grid*(knot[4]-knot[1]) - term4Grid*(knot[3]-knot[1]))/(knot[4]-knot[3])) / (knot[4]-knot[1])^2
u3Grid <- (term2Grid - (term3Grid*(knot[4]-knot[2]) - term4Grid*(knot[3]-knot[2]))/(knot[4]-knot[3])) / (knot[4]-knot[1])^2

true_y1_grid <- 100 + 1*u1Grid + 0*u2Grid + 0*u3Grid
true_y2_grid <- 100 + 0*u1Grid + 0*u2Grid + 1.5*u3Grid
true_y3_grid <- 100 + 1*u1Grid + 0*u2Grid + 1.5*u3Grid

# back-transform ests, SDs, and CI limits)
y2_Spline_Pred_Res$ests    <- y2_Spline_Pred_Res$ests    * 10 + 100
y2_Spline_Pred_Res$ciLower <- y2_Spline_Pred_Res$ciLower * 10 + 100
y2_Spline_Pred_Res$ciUpper <- y2_Spline_Pred_Res$ciUpper * 10 + 100
y2_Spline_Pred_Res$stderrs <- y2_Spline_Pred_Res$stderrs * 10

# evaluate performance
summarise_res(y2_Spline_Pred_Res, true_y2_grid)

# and with y3
y3_Spline_Pred_Res$ests    <- y3_Spline_Pred_Res$ests    * 10 + 100
y3_Spline_Pred_Res$ciLower <- y3_Spline_Pred_Res$ciLower * 10 + 100
y3_Spline_Pred_Res$ciUpper <- y3_Spline_Pred_Res$ciUpper * 10 + 100
y3_Spline_Pred_Res$stderrs <- y3_Spline_Pred_Res$stderrs * 10

summarise_res(y3_Spline_Pred_Res, true_y3_grid)
