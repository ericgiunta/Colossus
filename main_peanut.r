library(survival)
library(parallel)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppEigen)


Rcpp::sourceCpp("PEANUT_MODEL.cpp") 


print("RAN") #just to verify the R section has started
#----------------------------------------------------------------------------------------------------------_#
colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
df <- fread('../Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv',nThread=10,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)
df <- data.frame(sapply( df, as.numeric ))
df <- as.data.table(df)
setkey(df, age_exit, lung, age_entry)
#
stime <- Sys.time()
dfend <- df[lung==1, ]
#df <- df[age_exit>=min(dfend$age_exit),]
tu <- unique(dfend[,age_exit])
close( file( 'test.txt', open="w" ) ) #needs a clean file
print(length(tu))
for (t0 in tu){
    all_ind <- df[age_entry<t0 & age_exit>=t0, which = TRUE]
    event_ind <- df[age_exit==t0 & lung==1, which = TRUE]
    #
    cons_all <- c()
    for (i in all_ind){
        if (length(cons_all)==0){
            cons_all <- c(cons_all,i,i)
        } else if (cons_all[length(cons_all)]+1<i){
            cons_all <- c(cons_all,i,i)
        } else{
            cons_all <- c(cons_all[1:length(cons_all)-1],i)
        }
    }
    #
    cat(c(cons_all,event_ind[1],event_ind[length(event_ind)],'\n'),file='test.txt',sep=',',append=TRUE)
}
etime <- Sys.time()
print(etime-stime)
#
covariates <- c("cumulative_dose_lung_lag10","sexm","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
events <- c('lung')
ce <- c('lung','age_entry','age_exit')
#
a=rep(-.01,length(covariates))
a_lin=c(0.0,c())
a_loglin=c(0.0,a)
a_plin=c(0.0,c())
xe=as.matrix(df)
x_lin=as.matrix(df[,..events])
x_loglin=as.matrix(df[,..covariates])
x_plin=as.matrix(df[,..events])
term_bool=c(0,1,0)
theta=0.5
fir=1
modelform='M'
#
etime <- Sys.time()
print(etime-stime)
e <- peanut_transition(a_lin, a_loglin, a_plin,  x_lin,  x_loglin,  x_plin, fir, modelform,length(tu),term_bool)
print("all results")
print(e)
print(exp(100*e[1]))
res <- coxph(Surv(time=age_entry, time2=age_exit, event=lung) ~ cumulative_dose_lung_lag10 + sexm + YOB1 + YOB2 + YOB3 + YOB4 + COH_EDUC1 + COH_EDUC2 + COH_EDUC3 + COH_EDUC4 + COH_EDUC5 + COH_EDUC6 + COH_EDUC7 + COH_EDUC8 + COH_EDUC9, data =  df)
coefs <- unlist(res["coefficients"], use.names=FALSE)
print(coefs)
print(res["loglik"])
##
dff <- df[sexm==0,]
dfm <- df[sexm==1,]
##
stime <- Sys.time()
dfend <- dff[lung==1, ]
tu <- unique(dfend[,age_exit])
close( file( 'test.txt', open="w" ) ) #needs a clean file
print(length(tu))
for (t0 in tu){
    all_ind <- df[age_entry<t0 & age_exit>=t0, which = TRUE]
    event_ind <- df[age_exit==t0 & lung==1, which = TRUE]
    #
    cons_all <- c()
    for (i in all_ind){
        if (length(cons_all)==0){
            cons_all <- c(cons_all,i,i)
        } else if (cons_all[length(cons_all)]+1<i){
            cons_all <- c(cons_all,i,i)
        } else{
            cons_all <- c(cons_all[1:length(cons_all)-1],i)
        }
    }
    #
    cat(c(cons_all,event_ind[1],event_ind[length(event_ind)],'\n'),file='test.txt',sep=',',append=TRUE)
}
etime <- Sys.time()
print(etime-stime)
#
covariates <- c("cumulative_dose_lung_lag10","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
events <- c('lung')
ce <- c('lung','age_entry','age_exit')
#
a=rep(-.01,length(covariates))
a_lin=c(0.0,c())
a_loglin=c(0.0,a)
a_plin=c(0.0,c())
term_bool <- c(0,1,0)
xe=as.matrix(dff)
x_lin=as.matrix(dff[,..events])
x_loglin=as.matrix(dff[,..covariates])
x_plin=as.matrix(dff[,..events])
theta=0.5
fir=1
modelform='M'
#
etime <- Sys.time()
print(etime-stime)
e <- peanut_transition(a_lin, a_loglin, a_plin,  x_lin,  x_loglin,  x_plin, fir, modelform,length(tu),term_bool)
print("male results")
print(e)
print(exp(100*e[1]))
res <- coxph(Surv(time=age_entry, time2=age_exit, event=lung) ~ cumulative_dose_lung_lag10 + YOB1 + YOB2 + YOB3 + YOB4 + COH_EDUC1 + COH_EDUC2 + COH_EDUC3 + COH_EDUC4 + COH_EDUC5 + COH_EDUC6 + COH_EDUC7 + COH_EDUC8 + COH_EDUC9, data =  dff)
coefs <- unlist(res["coefficients"], use.names=FALSE)
print(coefs)
print(res["loglik"])
##
stime <- Sys.time()
dfend <- dfm[lung==1, ]
tu <- unique(dfend[,age_exit])
close( file( 'test.txt', open="w" ) ) #needs a clean file
print(length(tu))
for (t0 in tu){
    all_ind <- df[age_entry<t0 & age_exit>=t0, which = TRUE]
    event_ind <- df[age_exit==t0 & lung==1, which = TRUE]
    #
    cons_all <- c()
    for (i in all_ind){
        if (length(cons_all)==0){
            cons_all <- c(cons_all,i,i)
        } else if (cons_all[length(cons_all)]+1<i){
            cons_all <- c(cons_all,i,i)
        } else{
            cons_all <- c(cons_all[1:length(cons_all)-1],i)
        }
    }
    #
    cat(c(cons_all,event_ind[1],event_ind[length(event_ind)],'\n'),file='test.txt',sep=',',append=TRUE)
}
etime <- Sys.time()
print(etime-stime)
#
#covariates <- c("cumulative_dose_lung_lag10","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
#events <- c('lung')
#ce <- c('lung','age_entry','age_exit')
#
#a=rep(.01,length(covariates))
#a_lin=c(0.0,c())
#a_loglin=c(0.0,rep(.01,length(covariates)))
#a_plin=c(0.0,c())
xe=as.matrix(dfm)
x_lin=as.matrix(dfm[,..events])
x_loglin=as.matrix(dfm[,..covariates])
x_plin=as.matrix(dfm[,..events])
#theta=0.5
#fir=1
#modelform='M'
#
etime <- Sys.time()
print(etime-stime)
e <- peanut_transition(a_lin, a_loglin, a_plin,  x_lin,  x_loglin,  x_plin, fir, modelform,length(tu),term_bool)
print("female results")
print(e)
print(exp(100*e[1]))
res <- coxph(Surv(time=age_entry, time2=age_exit, event=lung) ~ cumulative_dose_lung_lag10 + YOB1 + YOB2 + YOB3 + YOB4 + COH_EDUC1 + COH_EDUC2 + COH_EDUC3 + COH_EDUC4 + COH_EDUC5 + COH_EDUC6 + COH_EDUC7 + COH_EDUC8 + COH_EDUC9, data =  dfm)
coefs <- unlist(res["coefficients"], use.names=FALSE)
print(coefs)
print(res["loglik"])
