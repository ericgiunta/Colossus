library(Colossus)
library(data.table)
library(parallel)


tfile <- file(paste(tempfile(), ".txt",sep="" ),open = "wt")
sink(file=tfile)
fname <- 'll_0.csv'
colTypes <- c( "double", "double", "double", "integer", "integer" )
df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
time1 <- "t0"
df$pyr <- df$t1-df$t0
pyr <- "pyr"
event <- "lung"
set.seed(3742)
df$rand <- floor(runif(nrow(df), min=0, max=5))
names <- c( "dose", "rand", "rand" )
term_n <- c(2,1,0)
tform <- c( "loglin", "lin", "plin" )
keep_constant <- c(0,0,0)
a_n <- c(0.01,0.1,0.1)
modelform <- "PAE"
fir <- 0
der_iden <- 0
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiters' = c(1,1), 'halfmax' = 5, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=4, 'ties'='breslow', 'double_step'=0)
strat_col <- "fac"

verbose <- FALSE
j_iterate <- 1
LL_comp <- c(-4.5574)
for (i in c(TRUE)){
    for (j in c(TRUE)){
        model_control <- list( 'strata'=i, 'single'=j)
        if (verbose){print(model_control)}
        a_n <- c(0.01,0.1,0.1)
        modelform <- "PAE"
        e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col,model_control)
#        expect_equal(e$LogLik,LL_comp[j_iterate],tolerance=1e-2)
    }
}
sink(NULL)
close(tfile)
