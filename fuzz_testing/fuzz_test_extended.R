# install.packages(c("fuzzr","Colossus","data.table"))

library(fuzzr)
library(Colossus)
library(data.table)


df <- data.table::data.table(
  UserID = c(112, 114, 213, 214, 115, 116, 117),
  Starting_Age = c(18, 20, 18, 19, 21, 20, 18),
  Ending_Age = c(30, 45, 57, 47, 36, 60, 55),
  Cancer_Status = c(0, 0, 1, 0, 1, 0, 0),
  a_val = c(0, 1, 1, 0, 1, 0, 1),
  b_val = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
  c_val = c(10, 11, 10, 11, 12, 9, 11),
  d_val = c(0, 0, 0, 1, 1, 1, 1),
  e_val = c(0, 0, 1, 0, 0, 0, 1)
)
control <- list(
  ncores = 1, lr = 0.75, maxiters = c(1, 1),
  halfmax = 1
)
formula_cox <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
  loglinear(a_val, b_val, c_val, 0) + plinear(d_val, 0) + multiplicative()

formula_casecon <- CaseCon_Strata(a_val, Cancer_Status) ~ loglinear(b_val, c_val) + plinear(d_val) + m()

formula_pois <- Pois(Ending_Age, Cancer_Status) ~ loglinear(CONST, a_val)

cox_res <- CoxRun(formula_cox, df, ncores = 1)
pois_res <- PoisRun(formula_pois, df, ncores = 1)

print(pois_res)
print(cox_res)
stop()

# test_all(), test_char(), test_int(), test_dbl(), test_lgl(), test_fctr(), test_date(), test_raw(), test_df(), test_null()



### ------------------- basic plot_option ------------------------ ###

fr <- fuzz_function(plot, "plot_options", x = cox_res, df = df, tests = test_all(), check_args = FALSE)
for (i in 1:length(fr)) {
  res <- fr[[i]]$test_result$errors
  test_name <- fr[[i]]$test_name
  test_name <- test_name[["x"]]
  if ((is.null(res)) && (!test_name %in% c())) {
    stop(paste("cox_res", test_name, sep = " "))
  }
}

### ------------------- individual plot options ------------------------ ###
# "verbose", "type", "age_unit", "strat_haz", "strat_col", "martingale", "km", "time_lims", "cov_cols", "studyid", "boundary"
bool_opts <- c("martingale", "km", "strat_haz")
str_opts <- c("type", "age_unit", "strat_col", "cov_cols", "studyid")
num_opts <- c("boundary", "verbose", "time_lims")

for (bool_arg in bool_opts) {
  fr <- fuzz_function(plot, bool_arg, x = cox_res, df = df, type = "risk", tests = test_all(), check_args = FALSE)
  for (i in 1:length(fr)) {
    res <- fr[[i]]$test_result$errors
    if (res == "missing value where TRUE/FALSE needed") {
        print(fr[[i]])
        stop()
    }
    test_name <- fr[[i]]$test_name
    test_name <- test_name[[bool_arg]]
    if ((is.null(res)) && (!test_name %in% c())) {
      print(fr[[i]])
      stop(paste(bool_arg, test_name, sep = " "))
    }
  }
}
