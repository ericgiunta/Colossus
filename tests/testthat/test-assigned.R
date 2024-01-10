## ------------------------------------- ##
## Verify working results
## ------------------------------------- ##
test_that("Poisson Assigned Events, no error", {
    df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
           "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
             "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
          "Cancer_Status"=c(12, 10, 18, 6, 1, 11, 4),
                      "a"=c(0,   1,   1,   0,   1,   0,   1),
                      "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
                      "c"=c(10,  11,  10,  11,  12,  9,   11),
                      "d"=c(0,   0,   0,   1,   1,   1,   1))

	df$pyr <- df$Ending_Age - df$Starting_Age
	pyr <- 'pyr'
	event <- "Cancer_Status"
	names <- c('a','b','c','d')
	Term_n <- c(0,1,1,2)
	tform <- c("loglin","lin","lin","plin")
	modelform <- "M"
	fir <- 0
	a_n <- c(-0.75, 0.1, -0.05, -1.5)

	keep_constant <- c(0,0,0,0)
	der_iden <- 0

	control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,'epsilon' = 1e-3,
	   'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
	   'dose_abs_max'=100.0,'verbose'=T, 'double_step'=1)
    #
    expect_no_error(RunPoissonEventAssignment(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that("Poisson Assigned Events, check results", {
    df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
           "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
             "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
          "Cancer_Status"=c(12, 10, 18, 6, 1, 11, 4),
                      "a"=c(0,   1,   1,   0,   1,   0,   1),
                      "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
                      "c"=c(10,  11,  10,  11,  12,  9,   11),
                      "d"=c(0,   0,   0,   1,   1,   1,   1))

	df$pyr <- df$Ending_Age - df$Starting_Age
	pyr <- 'pyr'
	event <- "Cancer_Status"
	names <- c('a','b','c','d')
	Term_n <- c(0,1,1,2)
	tform <- c("loglin","lin","lin","plin")
	modelform <- "M"
	fir <- 0
	a_n <- c(-0.75, 0.1, -0.05, -1.5)

	keep_constant <- c(0,0,0,0)
	der_iden <- 0

	control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,'epsilon' = 1e-3,
	   'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
	   'dose_abs_max'=100.0,'verbose'=FALSE, 'double_step'=1)
    #
    e <- RunPoissonEventAssignment(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
    expect_equal(e$caused,c(105.81165, -43.81165),tolerance=1)
    expect_equal(e$predict,c(134.79452, -53.34884),tolerance=1)
})
test_that("Cox Assigned Events, no error", {
    df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
		       "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
		         "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
		      "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
		                  "a"=c(0,   1,   1,   0,   1,   0,   1),
		                  "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
		                  "c"=c(10,  11,  10,  11,  12,  9,   11),
		                  "d"=c(0,   0,   0,   1,   1,   1,   1))
	# For the interval case
	time1 <- "Starting_Age"
	time2 <- "Ending_Age"
	event <- "Cancer_Status"
	names <- c('a','b','c','d')
	Term_n <- c(0,1,1,2)
	tform <- c("loglin","lin","lin","plin")
	modelform <- "M"
	fir <- 0
	a_n <- c(0.1, 0.1, 0.1, 0.1)

	keep_constant <- c(0,0,0,0)
	der_iden <- 0

	control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 50,'halfmax' = 5,
	   'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
	   'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
	   'verbose'=T, 'ties'='breslow','double_step'=1)

	expect_no_error(RunCoxEventAssignment(df, time1, time2, event, names, Term_n, tform,
		                 keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that("Poisson Assigned Events, check results", {
    df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
		       "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
		         "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
		      "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
		                  "a"=c(0,   1,   1,   0,   1,   0,   1),
		                  "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
		                  "c"=c(10,  11,  10,  11,  12,  9,   11),
		                  "d"=c(0,   0,   0,   1,   1,   1,   1))
	# For the interval case
	time1 <- "Starting_Age"
	time2 <- "Ending_Age"
	event <- "Cancer_Status"
	names <- c('a','b','c','d')
	Term_n <- c(0,1,1,2)
	tform <- c("loglin","lin","lin","plin")
	modelform <- "M"
	fir <- 0
	a_n <- c(0.1, 0.1, 0.1, 0.1)

	keep_constant <- c(0,0,0,0)
	der_iden <- 0

	control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 50,'halfmax' = 5,
	   'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
	   'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
	   'verbose'=F, 'ties'='breslow','double_step'=1)

	e <- RunCoxEventAssignment(df, time1, time2, event, names, Term_n, tform,
		                 keep_constant, a_n, modelform, fir, der_iden, control)
    expect_equal(e$caused,c(0.1447365, 1.8552635),tolerance=1e-2)
})
