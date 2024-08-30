library(Colossus)
library(data.table)
library(parallel)


# tfile <- file(paste(tempfile(), ".txt",sep="" ),open = "wt")
# sink(file=tfile)
fname <- 'll_comp_0.csv'
colTypes <- c( "double", "double", "double", "integer", "integer" )
df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
set.seed(3742)
df$rand <- floor(runif(nrow(df), min=0, max=5))

time1 <- "t0"
time2 <- "t1"


#
event <- "lung"
names <- c( "rand", "fac", "dose" )
term_n <- c(0,0,1)
tform <- c( "loglin", "loglin", "loglin" )
keep_constant <- c(0,0,0)
a_n <- c(-0.1,0.1,0.2)
modelform <- "M"
fir <- 0
der_iden <- 0

cens_weight <- c(0)

verbose <- FALSE

devs <- c()

modelform <- "M"
model_control <- list( 'strata'=FALSE, 'basic'=FALSE, 'single'=FALSE, 'cr'=FALSE)
a_n <- c(0.6465390, 0.4260961, 0.1572781)
lung_temp <- df$lung
df$lung <- rep(0,length(lung_temp))
keep_constant <- c(0,0,0)
print(df$lung)
#
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiters' = c(1,1), 'halfmax' = 2, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=4, 'ties'='breslow', 'double_step'=1)
RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="fac", model_control=model_control)
keep_constant <- c(1,1,1)
# sink(NULL)
# close(tfile)
