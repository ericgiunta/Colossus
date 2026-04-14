library(Colossus)
library(data.table)
df <- data.table::data.table(
    'UserID' = c(112, 114, 213, 214, 115, 116, 117),
    'Starting_Age' = c(18, 20, 18, 19, 21, 20, 18),
    'Ending_Age' = c(30, 45, 57, 47, 36, 60, 55),
    'Cancer_Status' = c(12, 10, 18, 6, 1, 11, 4),
    'a' = c(0, 1, 1, 0, 1, 0, 1),
    'b' = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
    'c' = c(10, 11, 10, 11, 12, 9, 11),
    'd' = c(0, 0, 0, 1, 1, 1, 1)
)

df$pyr <- df$Ending_Age - df$Starting_Age
a_n <- c(-0.75, 0.1, -0.05, -1.5)
keep_constant <- c(0, 0, 0, 0)
control <- list(
"ncores" = 1, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
"deriv_epsilon" = 1e-3, "step_max" = 0.2, "change_all" = TRUE,
"thres_step_max" = 100.0, "verbose" = 4
)
#
print("Poisson Run Start")
poisres <- PoisRun(Pois(pyr, Cancer_Status) ~ loglinear(a, 0) + linear(b, c, 1) + plinear(d, 2), df, a_n = a_n, control = control, norm = "max")
assign_control <- list(check_num = 4, verbose = 4)
print("Assignment Run Start")
e <- EventAssignment(poisres, df, assign_control = assign_control, z = 2)

elow <- e$lower_limit$predict
emid <- e$midpoint$predict
eupp <- e$upper_limit$predict
