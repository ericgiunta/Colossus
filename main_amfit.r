library(parallel)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppEigen)
library(parglm)

Rcpp::sourceCpp("AMFIT_MODEL.cpp") 




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
colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
df <- fread('../Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv',nThread=10,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)

df <- as.data.table(df) #verifys its a data.table
setkey(df, age_exit, lung, age_entry) #sorts by ages for ease of summarizing
df$Pyr <- df$age_exit -df$age_entry #defines person-year column
##
df <- within(df, {
    #list of new columns based on Linda's code
    lpy10k <- log((age_exit-age_entry)/10000)

    ## Sex indicators
    male <- sexm==0
    female <- sexm==1
    ##  msex <- 2*sex-1 # -1 for males; 1 for females

    ## Sex-specific age terms
    bc30 <- (YOB - 1930)/10 # birth cohort, centered at 1930, in decades
    lage70 <- log(age_entry/70) # Attained age, centered at 70 years
    lage70sq <- lage70^2 # Attained age squared
    lage70qsp <- lage70sq*(age_entry > 70) # Attained age quadratic spline with knot at 70 years
    mlage70 <- male*lage70
    flage70 <- female*lage70
    mlage70sq <- male*lage70sq
    flage70sq <- female*lage70sq
    mlage70qsp <- male*lage70qsp
    flage70qsp <- female*lage70qsp
    dosegy <- cumulative_dose_lung_lag10/1000  

    ## Sex-specific dose terms

    mdgy <- male*dosegy
    fdgy <- female*dosegy
})
# combine educational categories

df$educatcom <- ifelse (df$EDUC==1 | 
                 df$EDUC==3 | 
                 df$EDUC==4 | 
                 df$EDUC==6,1,NA)
df$educatcom <- ifelse (df$EDUC==2 | 
                 df$EDUC==5 | 
                 df$EDUC==7,2,df$educatcom)
df$educatcom <- ifelse (df$EDUC==8 | df$EDUC==9 | df$EDUC==10,3,df$educatcom)
## The original formula replicated
#        lung ~ -1 + factor(sexm) + factor(educatcom)
#          + bc30
#          + mlage70 + mlage70sq + mlage70qsp
#          + flage70 + flage70sq + flage70qsp
#          + offset(lpy10k)
#        df <- df[, lapply(.SD,sum), by=.(sexm,educatcom,bc30,mlage70 , mlage70sq , mlage70qsp,flage70 , flage70sq , flage70qsp), .SDcols=c("Pyr",'lung')]
#
data <- factorize(df,c('educatcom')) #splits into factors
edu_col <- data$cols
df <- data$df
#
covariates <- c("dosegy","sexm",edu_col,"bc30","mlage70" , "mlage70sq" , "mlage70qsp","flage70" , "flage70sq" , "flage70qsp") #complete list of covariates used
#
df$i_lin <- 1.0
df$i_loglin <- 1.0 #defines columns for intercept terms
df$i_plin <- 1.0
#
algorithm_control <- parglm.control(method = "FAST",trace=FALSE,nthreads=10) #control for the validated run
poisFit <- parglm(as.formula("lung ~ -1 + dosegy + sexm + educatcom + bc30 + mlage70 + mlage70sq + mlage70qsp + flage70 + flage70sq + flage70qsp"), family=poisson, data=df, control = algorithm_control) #contains the values to compare against
covariates_lin <- c('i_loglin')
covariates_loglin <- c("dosegy","sexm",edu_col,"bc30","mlage70" , "mlage70sq" , "mlage70qsp","flage70" , "flage70sq" , "flage70qsp") #Covariate list for each type of term
covariates_plin <- c('i_plin') #empty data.frames cause issues in the c++ section
ce <- c('Pyr','lung') #The duration and event columns
#
a_lin=c(0.0,c())
a_loglin=c(0.0,rep(-.01,length(covariates_loglin))) #the lists of covariates
a_plin=c(0.0,c()) #empty first value added to avoid empty vectors
term_bool <- c(0,1,0) #gives which terms are actually needed

x_lin=as.matrix(df[,..covariates_lin])
x_loglin=as.matrix(df[,..covariates_loglin]) #creates a matrix for each term type
x_plin=as.matrix(df[,..covariates_plin])
xe=as.matrix(df[,..ce]) #matrix of durations and events
fir=1 #some model types have a term that the rest are multiplied by, this gives the 0-indexed position of it. Required
modelform='M' #Model type, 'M', 'A', 'PA', 'PAE' for multiplicative, additive, product additive, and product additive excess
e <- amfit_transition(a_lin, a_loglin, a_plin,  x_lin,  x_loglin,  x_plin, fir, modelform,xe,term_bool) #performs run, returns the optimal parameters
print(covariates_loglin)
print(e)
coefs <- poisFit["coefficients"]
print(coefs)
res <- (unlist(poisFit["fitted.values"], use.names=FALSE) - unlist(poisFit["y"], use.names=FALSE))
print(sum(res^2)/nrow(df)) #the variance in error for comparison
