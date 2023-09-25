library(data.table)
library(parallel)
library(Colossus)

a <- c(0,1,2,3,4,5,6)
b <- c(1.1,2.5,3.0,4.5,5.1,6.8,7.2)
c <- c(0,1,0,0,0,1,0)
d <- c(3,4,5,6,7,8,9)
e <- c(0,1,0,1,0,0,1)
for (i in 1:10){
    a <- c(a,a)
    b <- c(b,b)
    c <- c(c,rep(0,length(c)))
    d <- c(d,d)
    e <- c(e,e)
}
df <- data.table("a"=a,"b"=b,"c"=c,"d"=d,"e"=e)
pyr <- "b"
event <- "c"
names <- c("d")
Term_n <- c(0)
tform <- c("loglin")
keep_constant <- c(0)
a_n <- c(-0.1)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
RunPoissonRegression_STRATA(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,c("e"))