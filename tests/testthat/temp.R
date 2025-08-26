library(Colossus)
library(data.table)
library(parallel)

fname <- "base_example.csv"
df <- fread(fname)
pyr <- "exit"
event <- "event"
names <- c("dose0", "dose1")
term_n <- c(0, 1)
tform <- c("loglin", "plin")
keep_constant <- c(0, 0)
a_n <- c(-2.917, 0.06526)
modelform <- "M"

#
model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = FALSE, "alpha" = 0.1)

control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
poisres <- PoisRun(Pois(exit, event) ~ loglinear(dose0, 0) + plinear(dose1, 0) + multiplicative(), df, a_n = a_n, control = control, keep_constant = keep_constant)
poisres$beta_0 <- c(-2.917, 0.06526)
# print(poisres)
# stop()
alpha_list <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
v_lower <- c(0.05418212, 0.04221923, 0.03135497, 0.02635009, 0.0100537, -0.0002002051, -0.009355479, -0.02009767, -0.02748292)
v_upper <- c(0.07591063, 0.0881787, 0.09953124, 0.10471923, 0.1221753, 0.1333447168, 0.143415518, 0.15550048, 0.16386229)
for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
  alpha <- alpha_list[alpha_i]
  a_n <- c(-2.917, 0.06526)
  model_control <- list("basic" = FALSE, "maxstep" = 20, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 1, "manual" = FALSE)
  curve_control <- list("alpha" = alpha, "para_number" = 1, "bisect" = TRUE, "step_size" = 0.5, "maxstep" = 20, "manual" = FALSE)
  e <- LikelihoodBound(poisres, df, curve_control, control = control)
  #    e <- PoissonCurveSolver(df, pyr, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control)
  a <- e$Parameter_Limits
  print(e)
}
