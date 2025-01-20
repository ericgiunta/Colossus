library(Colossus)
library(parallel)
library(data.table)

fname <- "ll_comp_0.csv"
colTypes <- c("double", "double", "double", "integer", "integer")
df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
set.seed(3742)
df$rand <- floor(runif(nrow(df), min = 0, max = 5))
df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
time1 <- "t0"
time2 <- "t1"
# df$censor <- (df$lung==0)
df$lung <- (df$lung > 0)
# event <- "censor"
names <- c("dose", "rand")
term_n <- c(0, 0)
tform <- c("loglin", "loglin")
realization_columns <- matrix(c("rand1","rand0"), nrow = 1)
realization_index <- c("rand")
keep_constant <- c(0, 0)
a_n <- c(0, 0)
modelform <- "M"
fir <- 0
der_iden <- 0
cens_weight <- c(0)
#
event <- "lung"
a_n <- c(-0.1, -0.1)
keep_constant <- c(0, 0)
verbose <- FALSE
model_control <- list("mcml" = TRUE)
a_n <- c(-0.1, -0.1)
#sink("out0.txt")
control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 100, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
e <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")
print(e)
#sink(NULL)

realization_columns <- matrix(c("rand0","rand1"), nrow = 1)
a_n <- c(-0.1, -0.1)
#sink("out1.txt")
e <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")
print(e)
#sink(NULL)