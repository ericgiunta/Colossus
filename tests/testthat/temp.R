library(Colossus)
library(data.table)
library(parallel)
library(survival)

fname <- "ll_comp_0.csv"
colTypes <- c("double", "double", "double", "integer", "integer")
df <- fread(fname, data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
set.seed(3742)
df$rand <- floor(runif(nrow(df), min = 0, max = 5))

time1 <- "t0"
time2 <- "t1"
df$censor <- (df$lung == 0)
event <- "censor"
names <- c("dose", "fac")
term_n <- c(0, 0)
tform <- c("loglin", "loglin")
keep_constant <- c(1, 0)
a_n <- c(0, 0)
modelform <- "M"
fir <- 0
der_iden <- 0
control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
plot_options <- list("name" = paste(tempfile(), "run", sep = ""), "verbose" = FALSE, "studyid" = "studyid", "age_unit" = "years")
dft <- GetCensWeight(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)
#
#
t_ref <- dft$t
surv_ref <- dft$surv
t_c <- df$t1
cens_weight <- approx(t_ref, surv_ref, t_c, rule = 2)$y
df$weighting <- cens_weight
#
event <- "lung"
a_n <- c(-0.1, -0.1)
keep_constant <- c(0, 0)

control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 4, "ties" = "breslow", "double_step" = 1)

verbose <- FALSE
j_iterate <- 1
LL_comp <- c(-81.12079, -81.12079, -77.97632, -77.97632, -68.60263, -68.67108, -75.34028, -75.3691, -81.12079, -81.12079, -77.97632, -77.97632, -68.60263, -68.67108, -75.34028, -75.3691, -122.8909, -122.8909, -119.9814, -119.9814, -109.6211, -109.6742, -117.0147, -117.0539, -122.8909, -122.8909, -119.9814, -119.9814, -109.6211, -109.6742, -117.0147, -117.0539)
for (i in c( FALSE)) {
    for (j in c( FALSE)) {
      for (k in c(TRUE)) {
        for (l in c( FALSE)) {
          model_control <- list("strata" = i, "basic" = j, "gradient" = k, "cr" = l)
          if (verbose) {
            print(model_control)
          }
          a_n <- c(-0.1, -0.1)
          control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 4, "ties" = "breslow", "double_step" = 1)
          e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
          print(e$beta_0)
          j_iterate <- j_iterate + 1
          a_n <- c(-0.1, -0.1)
          control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 4, "ties" = "efron", "double_step" = 0)
          e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
          print(e$beta_0)
        }
      }
    }
}
