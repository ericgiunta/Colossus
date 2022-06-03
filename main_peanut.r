library(survival)
library(parallel)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppEigen)


Rcpp::sourceCpp("PEANUT_MODEL.cpp") 

#' Defines Time Dependent Parameters
#' \code{time_factor} uses user provided bins and a list of columns to define interaction terms and update the data.table.
#' Technically could be used to define interaction terms for a binned column and a list of non-binned values
#'
#' @param df a data.table containing the columns of interest
#' @param time_bins an array of bins to apply to the time column to split by
#' @param col_list an array of column names that should have interaction terms defined
#' @param time_col the column name of the time column the data is binned by
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#'

time_factor <- function(df,time_bins,col_list,time_col){
    cols <- c()
    for (i in 1:length(col_list)){
        col <- col_list[i]
        for (i in 1:(length(time_bins)-1)){
            newcol <- c(paste(col,i,sep="_T"))
            if (sum(((df[,..time_col]>=time_bins[i])&(df[,..time_col]<time_bins[i+1])))>0){
                df[, newcol] <- df[,..col]*((df[,..time_col]>=time_bins[i])&(df[,..time_col]<time_bins[i+1]))
                cols <- c(cols, newcol)
                #
                newcol <- c(paste(i,"_T",sep=""))
                df[, newcol] <- 1*((df[,..time_col]>=time_bins[i])&(df[,..time_col]<time_bins[i+1]))
                cols <- c(cols, newcol)
            }
        }
    }
    print(df)
    list('df'=df, 'cols'=cols)
}

#' Splits a parameter into factors
#' \code{factorize} uses user provided list of columns to define new parameter for each unique value and update the data.table.
#' Not for interaction terms
#'
#' @param df a data.table containing the columns of interest
#' @param col_list an array of column names that should have factor terms defined
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#'

factorize <-function(df,col_list){
    cols <- c()
    for (i in 1:length(col_list)){
        col <- col_list[i]
        x <- sort(unlist(as.list(unique(df[,..col])),use.names=FALSE))
    #    print(x)
        for (i in x){
            newcol <- c(paste(col,i,sep="_"))
            if (sum(df[,..col]==i)>0){
                df[, newcol] <- 1*(df[,..col]==i)
                cols <- c(cols, newcol)
            }
        }
    }
#    print(df)
    list('df'=df, 'cols'=cols)
}




print("RAN") #just to verify the R section has started
#----------------------------------------------------------------------------------------------------------_#
colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
df <- fread('../Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv',nThread=10,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)

df <- as.data.table(df) #verifies its a data.table
setkey(df, age_exit, lung, age_entry) #sorts for risk group searching
#
stime <- Sys.time()
dfend <- df[lung==1, ] #for list of unique end times
#df <- df[age_exit>=min(dfend$age_exit),] #would remove data that ends before any risk group
#df <- df[age_entry<=max(dfend$age_exit),] #would remove data that ends after every risk group
tu <- unique(dfend[,age_exit])
close( file( 'test.txt', open="w" ) ) #needs a clean file
print(length(tu))
for (t0 in tu){
    all_ind <- df[age_entry<t0 & age_exit>=t0, which = TRUE] #finds row numbers for every row in the risk group
    event_ind <- df[age_exit==t0 & lung==1, which = TRUE] #finds row numbers for every event in the risk group, sequential due to the ordering
    #
    cons_all <- c() #stores the row numbers in pairs containing everything inbetween
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
    cat(c(cons_all,event_ind[1],event_ind[length(event_ind)],'\n'),file='test.txt',sep=',',append=TRUE) #writes to file to be loaded later on, non-constant length
}
etime <- Sys.time()
print(etime-stime)
#
covariates <- c("cumulative_dose_lung_lag10","sexm","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9") #covariates of interest
events <- c('lung') #event of interest
ce <- c('lung','age_entry','age_exit') #event and time columns
#
a=rep(-.01,length(covariates))
a_lin=c(0.0,c())
a_loglin=c(0.0,a) #the lists of covariates
a_plin=c(0.0,c()) #empty first value added to avoid empty vectors

x_lin=as.matrix(df[,..covariates_lin])
x_loglin=as.matrix(df[,..covariates_loglin]) #creates a matrix for each term type
x_plin=as.matrix(df[,..covariates_plin]))
term_bool=c(0,1,0)  #gives which terms are actually needed

fir=1 #some model types have a term that the rest are multiplied by, this gives the 0-indexed position of it. Required
modelform='M' #Model type, 'M', 'A', 'PA', 'PAE' for multiplicative, additive, product additive, and product additive excess
#
etime <- Sys.time()
print(etime-stime)
e <- peanut_transition(a_lin, a_loglin, a_plin,  x_lin,  x_loglin,  x_plin, fir, modelform,length(tu),term_bool) #performs run and return optimal
print("all results")
print(e) #the parameters
print(exp(100*e[1])) #the relative risk at 100 mGy estimated
res <- coxph(Surv(time=age_entry, time2=age_exit, event=lung) ~ cumulative_dose_lung_lag10 + sexm + YOB1 + YOB2 + YOB3 + YOB4 + COH_EDUC1 + COH_EDUC2 + COH_EDUC3 + COH_EDUC4 + COH_EDUC5 + COH_EDUC6 + COH_EDUC7 + COH_EDUC8 + COH_EDUC9, data =  df) #the value to compare against
coefs <- unlist(res["coefficients"], use.names=FALSE)
print(coefs)
print(res["loglik"])
#-------------------------------------------------------------------------------------------------------------------------------------#
#
#
#                           Past this point the same is done for the male and female specific groups
#
#
#-------------------------------------------------------------------------------------------------------------------------------------#
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
xe=as.matrix(dfm)
x_lin=as.matrix(dfm[,..events])
x_loglin=as.matrix(dfm[,..covariates])
x_plin=as.matrix(dfm[,..events])
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
