library(parallel)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppEigen)
library(parglm)

Rcpp::sourceCpp("AMFIT_MODEL.cpp") 

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
df <- data.frame(sapply( df, as.numeric ))
df <- as.data.table(df,key='age_entry')
#        df <- df[order(age_entry,age_exit)]
df$Pyr <- df$age_exit -df$age_entry
##
df <- within(df, {
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
##
#        lung ~ -1 + factor(sexm) + factor(educatcom)
#          + bc30
#          + mlage70 + mlage70sq + mlage70qsp
#          + flage70 + flage70sq + flage70qsp
#          + offset(lpy10k)
#        df <- df[, lapply(.SD,sum), by=.(sexm,educatcom,bc30,mlage70 , mlage70sq , mlage70qsp,flage70 , flage70sq , flage70qsp), .SDcols=c("Pyr",'lung')]
#
data <- factorize(df,c('educatcom'))
edu_col <- data$cols
df <- data$df
#
covariates <- c("dosegy","sexm",edu_col,"bc30","mlage70" , "mlage70sq" , "mlage70qsp","flage70" , "flage70sq" , "flage70qsp")
#
df$i_lin <- 1.0
df$i_loglin <- 1.0
df$i_plin <- 1.0
#
algorithm_control <- parglm.control(method = "FAST",trace=FALSE,nthreads=10)
poisFit <- parglm(as.formula("lung ~ -1 + dosegy + sexm + educatcom + bc30 + mlage70 + mlage70sq + mlage70qsp + flage70 + flage70sq + flage70qsp"), family=poisson, data=df, control = algorithm_control)
covariates_lin <- c('i_loglin')
covariates_loglin <- c("dosegy","sexm",edu_col,"bc30","mlage70" , "mlage70sq" , "mlage70qsp","flage70" , "flage70sq" , "flage70qsp")
covariates_plin <- c('i_plin')
ce <- c('Pyr','lung')
#
a_lin=c(0.0,c())
a_loglin=c(0.0,rep(-.01,length(covariates_loglin)))
a_plin=c(0.0,c())
term_bool <- c(0,1,0)
xe=as.matrix(df)
x_lin=as.matrix(df[,..covariates_lin])
x_loglin=as.matrix(df[,..covariates_loglin])
x_plin=as.matrix(df[,..covariates_plin])
xe=as.matrix(df[,..ce])
fir=1
modelform='M'
e <- amfit_transition(a_lin, a_loglin, a_plin,  x_lin,  x_loglin,  x_plin, fir, modelform,xe,term_bool)
print(covariates_loglin)
print(e)
coefs <- poisFit["coefficients"]
print(coefs)
res <- (unlist(poisFit["fitted.values"], use.names=FALSE) - unlist(poisFit["y"], use.names=FALSE))
print(sum(res^2)/nrow(df))
