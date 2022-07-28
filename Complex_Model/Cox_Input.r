library(survival)
library(parallel)
library(doParallel)
library(foreach)
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

#' Defines Interactions
#' \code{interact_them} uses user provided interactions define interaction terms and update the data.table.
#' assumes interaction is "+" or "*"
#' CURRENTLY DOES NOT APPLY ANTI_ALIASING
#'
#' @param df a data.table containing the columns of interest
#' @param interactions array of strings, each one is of form "term1?\*\?term2" for term1 interaction of type \*\ with term2, "?" dlimits
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#'

interact_them <- function(df,interactions,new_names){
    cols <- c()
    for (i in 1:length(interactions)){
        interac <- interactions[i]
        formula <- unlist(strsplit(interac,"?"),use.name=FALSE)
        if (length(formula)!=3){
            print(paste("Iteration:",interac,"has incorrect length of",length(formula),"but should be 3."))
            stop()
        }
        newcol <- paste(formula[1],formula[2],formula[3],sep="")
        if (new_names[i]!=''){
            newcol <- new_names[i]
        }
        col1 <- formula[1]
        col2 <- formula[3]
        if (paste(formula[1],"?",formula[2],"?",formula[3],sep="") %in% interactions[i+1:length(interactions)]){
            "its duped"
        } else if (paste(formula[3],"?",formula[2],"?",formula[1],sep="") %in% interactions[i+1:length(interactions)]){
            "the reverse is duped"
        } else {
            if (formula[2]=="+"){
                df[, newcol] <- df[,..col1] + df[,..col2]
                cols <- c(cols, newcol)
            } else if (formula[2]=="*"){
                df[, newcol] <- df[,..col1] * df[,..col2]
                cols <- c(cols, newcol)
            } else {
                print(paste("Incorrect operation of",formula[2]))
                stop()
            }
        }
    }
    print(df)
    list('df'=df, 'cols'=cols)
}

Likelihood_Ratio_Test <- function(alternative_model, null_model){
    return -2*(alternative_model["LogLik"] - null_model["LogLik"])
}

DoseForm <- function(df, parameters){
    names = unlist(parameters['names'],use.name=FALSE)
    types = unlist(parameters['types'],use.name=FALSE)
    intercepts = unlist(parameters['intercepts'],use.name=FALSE)
    steps = unlist(parameters['steps'],use.name=FALSE)
    i=0
    cols<-c()
    d_types<-c()
#    print(names(df))
    for (i in 1:length(names)){
        type <- types[i]
        name <- names[i]
        if (type=="LNT"){
            new_name <- paste(name,"_l",i,sep="")
            df[,new_name] <- df[,..name]
            d_types <- c(d_types,"lin")
        } else if (type=="quadratic"){
            new_name <- paste(name,"_q",i,sep="")
            df[,new_name] <- df[,..name] * df[,..name]
            d_types <- c(d_types,"lin")
        } else if (type=="lin-quad"){
            print(paste("Linear-Quadratic isn't implemented"))
            stop()
        } else if (type=="lin-exp"){
            new_name <- paste(name,"_le",i,sep="")
            df[,new_name] <- df[,..name]
            d_types <- c(d_types,"loglin")
        } else if (type=="LT"){
            new_name <- paste(name,"_lt",i,sep="")
            df[,new_name] <- df[,..name] - intercepts[i]
            df[get(new_name)<0.0] <- 0.0
            d_types <- c(d_types,"lin")
        } else if (type=="step"){
            new_name <- paste(name,"_s",i,sep="")
            df[,new_name] <- 1
            df[get(new_name)<steps[i]] <- 0.0
            d_types <- c(d_types,"lin")
        } else if (type=="inv-step"){
            new_name <- paste(name,"_is",i,sep="")
            df[,new_name] <- 1
            df[get(new_name)>steps[i]] <- 0.0
            d_types <- c(d_types,"lin")
        } else if (type=="step-slope"){
            new_name <- paste(name,"_ss",i,sep="")
            df[,new_name] <- df[,..name] - intercepts[i]
            df[get(new_name)<steps[i]] <- 0.0
            d_types <- c(d_types,"lin")
        }
        cols <- c(cols,new_name)
        i=i+1
    }
    return (list('df'=df, 'names'=cols, 'd_types'=d_types))
#    return list('df'=df, 'names'=cols)
}


run_coxph <-function(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event){
    #-------------------------------------------------------------------------------------------------------------#
    #   df is the data that will be used
    #   df is changed to a data.table to make filtering easier
    colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)
    setkeyv(df, c(time2, event,time1))
    #-------------------------------------------------------------------------------------------------------------#
    #   The goal is to precompute the indices needed for each event time
    #   The file has structure {start of all data, end of all data, first event index, second event index, ...,}
    #
    #
    dfend <- df[get(event)==1, ]
    #
    tu <- unlist(unique(dfend[,..time2]))
    if (!unlist(control["run"],use.name=FALSE)){
        close( file( 'test.txt', open="w" ) ) #needs a clean file
        ce <- c(time1,time2,event)
        Write_Ind_File(as.matrix(df[,..ce]),tu)
    }
    #
    dose_n <- dose_paras$names
    #
    all_names <- c()
    all_names <- c(all_names, dose_n)
    all_names <- c(all_names, lin_n)
    all_names <- c(all_names, loglin_n)
    all_names <- c(all_names, plin_n)
#    print(df,nrows=10)
    #
    #-------------------------------------------------------------------------------------------------------------#
#    print(lin_n)
    if (length(lin_n)==0){
        lin_n = c(event)
    }
#    print(lin_n)
    if (length(loglin_n)==0){
        loglin_n = c(event)
    }
    if (length(plin_n)==0){
        plin_n = c(event)
    }
    x_lin=as.matrix(df[,..lin_n])
    x_loglin=as.matrix(df[,..loglin_n])
    x_plin=as.matrix(df[,..plin_n])
    x_dose=as.matrix(df[,..dose_n])
    term_bool=c(0,0,0)
    if (length(a_lin)>0){
        term_bool[1]=1
    }
    if (length(a_loglin)>0){
        term_bool[2]=1
    }
    if (length(a_plin)>0){
        term_bool[3]=1
    }
    ##
#    print(length(tu))
    print("all results")
    e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras)
    print(e)
    b = e$beta_0
    er = e$Standard_Error
    for (i in 1:length(b)){
        itemp = 1.645 * sqrt(er[i])
        temp = paste(all_names[i],": ",sprintf("%.5e",b[i])," +/- ",sprintf("%.5e",itemp)," (",sprintf("%.5e",b[i]-itemp),", ",sprintf("%.5e",b[i]+itemp),")",sep="")
        print(temp)
    }
}

#fname <- '../Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv'
fname <- '/home/user/Documents/CMS/Combined_MW_NPP_IR_Lung_Lag10_-_TEST_4.18.22.csv'
lin_n <- c()
loglin_n <- c("sexm","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
plin_n <- c()

a_lin=c()
a_loglin <- rep(-.01,length(loglin_n))
a_plin=c()

modelform <- 'M'
fir=0

control <- list('lr' = 0.5,'maxiter' = -1,'halfmax' = -1,'epsilon' = 1e-7,'dbeta_max' = 0.2,'deriv_epsilon' = 1e-5,'run'=TRUE,'batch_size'=10000000, 'abs_max'=1.0)

time1="age_entry"
time2="age_exit"
event="lung"

dose_paras <- list('names' =c('cumulative_dose_lung_lag10'), 'beta_loglin_slope'=list(c(1.0)), 'beta_loglin_top'=list(c(-.1)), 'beta_lin_slope'=list(c(.1)), 'beta_lin_int'=list(c(0.0)),'beta_quad'=list(c(.1)),'beta_step_slope'=list(c(0.0)),'beta_step_int'=list(c(0.0)))

#
#//df_dose :matrix of doses
#//beta_loglin_slope: list of loglin linear magnitudes
#//beta_loglin_top: list of loglin exponential magnitudes
#//beta_lin_slope: list of linear slopes
#//beta_lin_int: list of linear intercepts
#//beta_quad: list of quadratic slopes
#//beta_step_slope: list of step magnitudes
#//beta_step_int: list of step intercepts
#
beta_loglin_slope <- list(c(1.0))
beta_loglin_top <- list(c(-.1))
beta_lin_slope <- list(c(.1))
beta_lin_int <- list(c(0.0))
beta_quad <- list(c(.1))
beta_step_slope <- list(c(0.0))
beta_step_int <- list(c(0.0))
#
run_coxph(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event)


