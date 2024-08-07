library(data.table)
library(parallel)
library(Colossus)

fname <- 'dose.csv'
colTypes <- c( "double", "double", "double", "integer" )
df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
time1 <- "t0"
time2 <- "t1"
event <- "lung"
names <- c( "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose" )
term_n <- c(0,0,0,0,0,0,0,0,0,0,0)
tform <- c( "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope" )
keep_constant <- c(0,0,0,0,0,0,0,0,0,0,0)
tform <- c( "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope" )
a_n <-   c(-0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1)
modelform <- "M"
fir <- 0
der_iden <- 0
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = 20, 'halfmax' = 5, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
e <- RunCoxRegression(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
print(e)
stop()
fname <- 'll_0.csv'
colTypes <- c( "double", "double", "double", "integer", "integer" )
df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
time1 <- "t0"
df$pyr <- df$t1-df$t0
pyr <- "pyr"
event <- "lung"
set.seed(3742)
df$rand <- floor(runif(nrow(df), min=0, max=5))
names <- c( "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose",  "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand" )
term_n <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1)
tform <- c( "loglin_slope", "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope", "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope" )
keep_constant <- c(0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
a_n <-   c(1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1)

modelform <- "PAE"
fir <- 0
der_iden <- 0
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiters' = c(1,1), 'halfmax' = 5, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=0)
strat_col <- "fac"

verbose <- FALSE
j_iterate <- 1
LL_comp <- c(-496.736611833105, -475.421339181129, -466.96134309201, -463.218078300204, -4497.1782941901, -3577.95315136081, -2557.93228186619, -2350.19459732101)
for (i in c(TRUE,FALSE)){
    for (j in c(TRUE,FALSE)){
        model_control <- list( 'strata'=i, 'single'=j)
        if (verbose){print(model_control)}
        a_n <-   c(1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1)
        modelform <- "PAE"
        e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col,model_control)
        print(paste(e$LogLik,LL_comp[j_iterate],sep=", "))
        j_iterate <- j_iterate + 1
        a_n <-   c(1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1)
        modelform <- "A"
        e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col,model_control)
        print(paste(e$LogLik,LL_comp[j_iterate],sep=", "))
        j_iterate <- j_iterate + 1
        if (verbose){print( "---------------" )}
    }
}
