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
    df <- fread(fname,nThread=detectCores(),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)
    setkeyv(df, c(time2, event,time1))
    #-------------------------------------------------------------------------------------------------------------#
    #   The goal is to precompute the indices needed for each event time
    #   The file has structure {start of all data, end of all data, first event index, second event index, ...,}
    dfend <- df[get(event)==1, ]
    #
    tu <- unlist(unique(dfend[,..time2]))
    if (!unlist(control["run"],use.name=FALSE)){
        close( file( 'test.txt', open="w" ) ) #needs a clean file
        for (t0 in tu){
            all_ind <- df[(get(time1)<t0 & get(time2)>=t0), which = TRUE] #finds row numbers for every row in the risk group
            event_ind <- df[(get(time2)==t0 & get(event)==1), which = TRUE] #finds row numbers for every event in the risk group, sequential due to the ordering
#            print(df[(get(time1)<t0 & get(time2)>=t0), which = TRUE])
#            print(df[(get(time2)==t0 & get(event)==1), which = TRUE])
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
    }
    #
#    print(names(df))
    d_control <- DoseForm(df, dose_paras)
    doseform <- dose_paras$doseform
    df <- d_control$df
    dose_n <- d_control$names
    dose_terms <- d_control$d_types
    #
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
    lr = control$lr
    maxiter = control$maxiter
    halfmax = control$halfmax
    epsilon = control$epsilon
    dbeta_max = control$dbeta_max
    deriv_epsilon = control$deriv_epsilon
    batch_size = control$batch_size
    ##
    for (i in 1:2000){
        print(a_{term})
        e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),a_dose,  x_lin,  x_loglin,  x_plin,x_dose, fir, modelform,doseform,dose_terms,length(tu),term_bool, lr, maxiter, halfmax, epsilon, dbeta_max, deriv_epsilon,batch_size)
        a_{term}[{ind2}] <- runif(1,{spread}) + {best}
}

#fname <- '../Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv'
fname <- '../Documents/CMS/Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv'
{term1}_n <- c() <- c()
{term}_n <- c("sexm","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
{term0}_n <- c()

a_{term1}=c()
a_{term} <- c(-0.0006300516,-0.0630181076,-0.3932366917,-0.6613753288,-1.1079753794,-0.7946066410,-0.1990600045,0.1616788011,-0.6317184838,-0.0572978471,-0.5560610214,-1.1076729774,-1.2444050074,-1.4178535299)
a_{term0}=c()
a_dose=c(0.0001946145)

modelform <- '{model}'
fir={ind1}


control <- list('lr' = 0.95,'maxiter' = 1,'halfmax' = 1,'epsilon' = 1e-8,'dbeta_max' = 0.1,'deriv_epsilon' = 1e-8,'run'=TRUE,'batch_size'=1000)

time1="age_entry"
time2="age_exit"
event="lung"

dose_paras <- list('names' =c('cumulative_dose_lung_lag10'), 'types'=c('lin-exp'), 'intercepts'=c(0), 'steps'=c(0), 'doseform'='A')

run_coxph(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event)

