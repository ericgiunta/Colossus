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


get_conf_int <-function(alpha){
    q1 <- qchisq(1-alpha, df=1)

}

Check_Dupe_Columns <- function(df,cols){
    ##
    if (length(cols)>1){
        features_pair <- combn(cols, 2, simplify = F) # list all column pairs
        toRemove <- c() # init a vector to store duplicates
        for(pair in features_pair) { # put the pairs for testing into temp objects
          f1 <- pair[1]
          f2 <- pair[2]

          if (!(f1 %in% toRemove) & !(f2 %in% toRemove)) {
            if (all(df[[f1]] == df[[f2]])) { # test for duplicates
              print(paste(f1, " and ", f2, " are equals.
",sep=""))
              toRemove <- c(toRemove, f2) # build the list of duplicates
            }
          }
        }
        return(setdiff(cols, toRemove))
    } else if (length(cols)==1){
        if (min(df[,..cols])==max(df[,..cols])){
            return(c())
        } else {
            return(cols)
        }
    } else {
        return(c())
    }
    return(c())
}

Check_Trunc <- function(df,ce){
    if (ce[1]=="%trunc%"){
        tname <- ce[2]
        tmin <- min(df[,get(tname)])
        df[,':='(right_trunc=tmin)]
        ce[1]="right_trunc"
    } else if (ce[2]=="%trunc%") {
        tname <- ce[1]
        tmax <- max(df[,get(tname)])
        df[,':='(left_trunc=tmax)]
        ce[2]="left_trunc"
    }
    return (list('df'=df,'ce'=ce))
}


control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 1,'epsilon' = 1e-5,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=FALSE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='efron')

time1="%trunc%"
time2="time"
event="cases"

#
df <- read.table("LUNG.DAT",header=FALSE,skip=0,col.names=c("cases","time","status","dxmnths","age","prior","trmnt","cell"))
df <- as.data.table(df)
setkeyv(df, c(time2, event))

df[,stat50 := status - 50]
df[,age60 :=  age - 60]
df[,trmnt := trmnt - 1 ]
df[,prior := prior / 10 ]
df0 <- data.table("cell_0"=c(0), "cell_1"=c(0), "cell_2"=c(0), "cell_3"=c(0), "stat50"=c(40), "age60"=c(20), "dxmnths"=c(6), "prior"=c(1), "trmnt"=c(1)) 


val <- factorize(df, c("cell"))

df <- val$df

names <- c(val$cols,"stat50", "age60", "dxmnths", "prior", "trmnt")

#
modelform <- 'A'
fir=0
der_iden=3
Term_n <- rep(0,length(names))
#STerm_n <- c(1,0,3,2,4,6,5)
#STerm_n <- c(STerm_n,rep(0,length(Term_n)-length(STerm_n)))
tform <- rep("loglin",length(names))
a_n <- rep(-0.1,length(names))
keep_constant <- rep(0,length(names))

#-------------------------------------------------------------------------------------------------------------#
#   The goal is to precompute the indices needed for each event time
#   The file has structure {start of all data, end of all data, first event index, second event index, ...,}
#
#
dfend <- df[get(event)==1, ]
#
tu <- unlist(unique(dfend[,..time2]))
print(length(tu))
#
dfc <- c()
dmax <- 0
all_names <- unique(names)
dfc <- match(names,all_names)
#-------------------------------------------------------------------------------------------------------------#
##
term_tot <- max(Term_n)+1
x_all=as.matrix(df[,..all_names])
x0_all=as.matrix(df0[,..all_names])
##
#    print(length(tu))
ce <- c(time1,time2,event)
#
t_check <- Check_Trunc(df,ce)
df <- t_check$df
ce <- t_check$ce
#
print(df)
print(names)
print(ce)
term_tot <- max(Term_n)+1
out <- tryCatch(
    {
        #R <- peanut_risk_sub(Term_n, tform, a_n,dfc,x0_all, fir, modelform, control, term_tot)
        e <- peanut_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,..ce]),tu,keep_constant,term_tot)
        #
    },
    error=function(cond){message(cond)},
    warning=function(cond){message(cond)},
    finally={print("_________________")}
)
