library(survival)
library(parallel)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppEigen)


Rcpp::sourceCpp("../PEANUT_MODEL.cpp")


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
              print(paste(f1, " and ", f2, " are equals.\n",sep=""))
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
calc_coxph <-function(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event){
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
    print(length(tu))
    #
    if (length(lin_n)!=0){
        lin_n = Check_Dupe_Columns(df,lin_n)
    }
#    print(lin_n)
    if (length(loglin_n)!=0){
        loglin_n = Check_Dupe_Columns(df,loglin_n)
    }
    if (length(plin_n)!=0){
        plin_n = Check_Dupe_Columns(df,plin_n)
    }
    if (length(dose_paras$names)!=0){
        dose_paras$names = Check_Dupe_Columns(df,dose_paras$names)
    }
    #
    dose_n <- dose_paras$names
    #
    all_names <- c()
    all_names <- c(all_names, dose_paras$terms)
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
    ce <- c(time1,time2,event)
    print("all results")
    print(nrow(df))
    print(length(tu))
    e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu)
    print(e)
    b = e$beta_0
    er = e$Standard_Deviation
    for (i in 1:length(b)){
        temp = paste(all_names[i],": ",sprintf("%.5e",b[i])," +/- ",sprintf("%.5e",2.0*er[i])," (",sprintf("%.5e",b[i]-2.0*er[i]),", ",sprintf("%.5e",b[i]+2.0*er[i]),")",sep="")
        print(temp)
    }
    #
    loglin_n <- c("YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
    a_loglin <- rep(-.01,length(loglin_n))
    all_names <- c()
    all_names <- c(all_names, dose_paras$terms)
    all_names <- c(all_names, lin_n)
    all_names <- c(all_names, loglin_n)
    all_names <- c(all_names, plin_n)
    dfm <- df[sexm==0,]
    dff <- df[sexm==1,]
    dfend <- dfm[get(event)==1, ]
    #
    tu <- unlist(unique(dfend[,..time2]))
    print(length(tu))
    x_lin=as.matrix(dfm[,..lin_n])
    x_loglin=as.matrix(dfm[,..loglin_n])
    x_plin=as.matrix(dfm[,..plin_n])
    x_dose=as.matrix(dfm[,..dose_n])
    print("male results")
    print(nrow(dfm))
    print(length(tu))
    e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(dfm[,..ce]),tu)
    print(e)
    b = e$beta_0
    er = e$Standard_Deviation
    for (i in 1:length(b)){
        temp = paste(all_names[i],": ",sprintf("%.5e",b[i])," +/- ",sprintf("%.5e",2.0*er[i])," (",sprintf("%.5e",b[i]-2.0*er[i]),", ",sprintf("%.5e",b[i]+2.0*er[i]),")",sep="")
        print(temp)
    }
    dfend <- dff[get(event)==1, ]
    #
    tu <- unlist(unique(dfend[,..time2]))
    print(length(tu))
    x_lin=as.matrix(dff[,..lin_n])
    x_loglin=as.matrix(dff[,..loglin_n])
    x_plin=as.matrix(dff[,..plin_n])
    x_dose=as.matrix(dff[,..dose_n])
    print("female results")
    print(nrow(dff))
    print(length(tu))
    e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(dff[,..ce]),tu)
    print(e)
    b = e$beta_0
    er = e$Standard_Deviation
    for (i in 1:length(b)){
        temp = paste(all_names[i],": ",sprintf("%.5e",b[i])," +/- ",sprintf("%.5e",2.0*er[i])," (",sprintf("%.5e",b[i]-2.0*er[i]),", ",sprintf("%.5e",b[i]+2.0*er[i]),")",sep="")
        print(temp)
    }
}

plot_coxph <-function(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event){
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
    print(length(tu))
    print("its going")
    dose_n <- dose_paras$names
    #
    all_names <- c()
    all_names <- c(all_names, dose_paras$names)
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
    ce <- c(time1,time2,event)
    e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu)
    print(e)
    b = e$beta_0
    er = e$Standard_Deviation
#    print(er)
#    print("all results")
#    print(nrow(df))
#    print(length(tu))
#    print(names(dose_paras))
    Plot_Type=c("SURV","not_used")
    print("start survival calculations")
    #                      a_lin,        a_loglin,        a_plin,   x_lin,  x_loglin,  x_plin, x_dose, fir, modelform, ntime,    include_bool,  Control,  Dose_paras,  d  f_groups,  tu 
    e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,0)
    t <- c()
    ch <- c()
    surv <- c()
    dt <- 1
    dft=data.table("time"=tu,"base"=e$baseline)
    for (i in tu[1]:tu[length(tu)]){
        t <- c(t,i)
        temp <- sum(dft[time<i, base])
        ch <- c(ch, temp)
        surv <- c(surv, exp(-1*temp))
    }
    print("martin plots")
#    print(unique(df[,studyID]))
    time_s <- df[,get(time1)]
    time_e <- df[,get(time2)]
    ch_fun <- approxfun(x=t,y=ch,rule=2)
    ch_e <- ch_fun(time_e)
    ch_s <- ch_fun(time_s)
    #
#    stop()
    e_i <- df[,get(event)]
    dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"IDS"=df$studyID,"time"=time_e,"doses"=df$cumulative_dose_lung_lag10)
    dfr$res = dfr$e - (dfr$ch_e-dfr$ch_s) * dfr$Risks
    #
#    temp_list <- c("time","res")
    Martingale_Error <- dfr[, lapply(.SD,sum), by=IDS]
    times <- dfr[, lapply(.SD,max), by=IDS]
    ##
    print("plotting dose by dose")
    jpeg('dose_sum_by_max.jpg', units="in", width=5, height=5, res=1200)
#    smoothScatter(times[,doses],Martingale_Error[,doses], xlab="Max Dose",ylab="Sum of Doses",nbin=100, colramp= colorRampPalette(c("white","red")))
    plot(times[,doses],Martingale_Error[,doses], xlab="Max Dose",ylab="Sum of Doses",type='p')
    dev.off()
    ##
    print("plotting martingale residuals")
    jpeg('martin_plot.jpg', units="in", width=5, height=5, res=1200)
    smoothScatter(times[,time],Martingale_Error[,res], xlab="age",ylab="Martingale Residuals",nbin=100, colramp= colorRampPalette(c("white","red")))
#    plot(times[,time],Martingale_Error[,res], xlab="age",ylab="Martingale Residuals")
    dev.off()
    print("plotting martingale residuals by dose")
    jpeg('martin_plot_dose.jpg', units="in", width=5, height=5, res=1200)
    smoothScatter(times[doses>0,doses],Martingale_Error[doses>0,res], xlab="Cumulative Dose",ylab="Martingale Residuals",nbin=200, colramp= colorRampPalette(c("white","red")))
#    plot(times[,doses],Martingale_Error[,res], xlab="age",ylab="Martingale Residuals")
    dev.off()
#    #    stop()
#    #
    print("plotting cumulative hazard")
    jpeg('ch_plot.jpg', units="in", width=5, height=5, res=1200)
    plot(t,ch, xlab="age",ylab="Cumulative Hazard")
    dev.off()
    #
    print("plotting survival curve")
    jpeg('surv_plot.jpg', units="in", width=5, height=5, res=1200)
    plot(t,surv, type="s", xlab="age",ylab="Survival")
    dev.off()
    #
    fir=0
    for (fir in 1:length(all_names)){
        print(paste("survival:",all_names[fir]))
        lfir <- c(all_names[fir])
        dfend <- df[get(event)==1, ]
        fir_c <- mean(dfend[,get(all_names[fir])])
        fmean <- mean(fir_c)
        print(fmean)
        df_u <- df[get(all_names[fir])>=fmean,]
        df_l <- df[get(all_names[fir])<=fmean,]
        #
        u_num = unlist(unique(df_u[,studyID]),use.names=FALSE)
        l_num = unlist(unique(df_l[,studyID]),use.names=FALSE)
        #
        t_u <- c(0)
        t_l <- c(0)
        n_u <- c(1)
        n_l <- c(1)
        tu <- unlist(unique(dfend[,..time2]))
        for (i in tu[1]:tu[length(tu)]){
            #
            df0 <- df_u[get(time2)<i,]
            u_num = length(unlist(unique(df0[,studyID]),use.names=FALSE))
            u_ev = sum(df0[, get(event)])
            df0 <- df_l[get(time2)<i,]
            l_num = length(unlist(unique(df0[,studyID]),use.names=FALSE))
            l_ev = sum(df0[, get(event)])
            #
            if (u_num>0){
                t_u <- c(t_u,i)
                temp <- (u_num - u_ev)/u_num
                n_u <- c(n_u, temp)
            }
            if (l_num>0){
                t_l <- c(t_l,i)
                temp <- (l_num - l_ev)/l_num
                n_l <- c(n_l, temp)
            }
        }
        print(min(n_u))
        print(min(n_l))
        print("plotting KM")
        jpeg(paste("KM_",fir,".jpg",sep=""), units="in", width=5, height=5, res=1200)
        plot(t_u,n_u, xlab=all_names[fir],ylab="Survival Fraction",ylim=c(min(n_l),max(n_u)),col='red',type='l')
        lines(t_l,n_l,col='blue')
        legend("topright", legend=c("Below Mean", "Above Mean"), col=c("red", "blue"), lty=1:2, cex=0.8,title=all_names[fir])
        dev.off()
        x_lin=as.matrix(df_l[,..lin_n])
        x_loglin=as.matrix(df_l[,..loglin_n])
        x_plin=as.matrix(df_l[,..plin_n])
        x_dose=as.matrix(df_l[,..dose_n])
        print(dim(x_dose))
        dfend <- df_l[get(event)==1, ]
        #
        tu <- unlist(unique(dfend[,..time2]))
        e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df_l[,..ce]),tu,Plot_Type,0)
        t <- c()
        ch <- c()
        surv <- c()
        dt <- 1
        dft=data.table("time"=tu,"base"=e$baseline)
        for (i in tu[1]:tu[length(tu)]){
            t <- c(t,i)
            temp <- sum(dft[time<i, base])
            ch <- c(ch, temp)
            surv <- c(surv, exp(-1*temp))
        }
        #
        Ls <- log(surv)
        LLs_l <- log(-Ls)
        Lt_l <- log(t)
        ##
        x_lin=as.matrix(df_u[,..lin_n])
        x_loglin=as.matrix(df_u[,..loglin_n])
        x_plin=as.matrix(df_u[,..plin_n])
        x_dose=as.matrix(df_u[,..dose_n])
        dfend <- df_u[get(event)==1, ]
        #
        tu <- unlist(unique(dfend[,..time2]))
        e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df_u[,..ce]),tu,Plot_Type,0)
        t <- c()
        ch <- c()
        surv <- c()
        dt <- 1
        dft=data.table("time"=tu,"base"=e$baseline)
        for (i in tu[1]:tu[length(tu)]){
            t <- c(t,i)
            temp <- sum(dft[time<i, base])
            ch <- c(ch, temp)
            surv <- c(surv, exp(-1*temp))
        }
        #
        Ls <- log(surv)
        LLs_u <- log(-Ls)
        Lt_u <- log(t)
        ##
        jpeg(paste('log_log_surv_plot',fir,'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
        plot(Lt_l,LLs_l, xlab="Log-Age",ylab="Log of Log Survival",col='red',type='l')
        lines(Lt_u,LLs_u,col='blue')
        legend("topleft", legend=c("Below Mean", "Above Mean"), col=c("red", "blue"), lty=1:2, cex=0.8,title=all_names[fir])
        lines(Ls ~ t)
        dev.off()
    }
#    stop()
    Plot_Type=c("RISK","not_used")
    fir=0
    dfend <- df[get(event)==1, ]
    #
    tu <- unlist(unique(dfend[,..time2]))
    print("start risk calculations")
    for (fir in 1:length(all_names)){
        print(all_names[fir])
        #
        dose_paras <- list('names' =c('cumulative_dose_lung_lag10'),'terms'=c('beta_loglin_top'), 'beta_loglin_slope'=list(c(1.0)), 'beta_loglin_top'=list(c(0.0001946145)), 'beta_lin_slope'=list(c(0.0)), 'beta_lin_int'=list(c(0.0)),'beta_quad'=list(c(0.0)),'beta_step_slope'=list(c(0.0)),'beta_step_int'=list(c(0.0)))
        a_loglin <- c(-0.0006300516,-0.0630181076,-0.3932366917,-0.6613753288,-1.1079753794,-0.7946066410,-0.1990600045,0.1616788011,-0.6317184838,-0.0572978471,-0.5560610214,-1.1076729774,-1.2444050074,-1.4178535299)
        if (fir==1){
            #
            lfir <- c(all_names[fir])
            uniq <- unlist(unique(df[,..lfir]))
    #        print(uniq)
            #
            e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir-1, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,length(uniq))
            x <- e$x
            y <- e$y
            #
            dose_paras <- list('names' =c('cumulative_dose_lung_lag10'),'terms'=c('beta_loglin_top'), 'beta_loglin_slope'=list(c(1.0)), 'beta_loglin_top'=list(b[1]-er[1]), 'beta_lin_slope'=list(c(0.0)), 'beta_lin_int'=list(c(0.0)),'beta_quad'=list(c(0.0)),'beta_step_slope'=list(c(0.0)),'beta_step_int'=list(c(0.0)))
            e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir-1, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,length(uniq))
            xl <- e$x
            yl <- e$y
            #
            dose_paras <- list('names' =c('cumulative_dose_lung_lag10'),'terms'=c('beta_loglin_top'), 'beta_loglin_slope'=list(c(1.0)), 'beta_loglin_top'=list(b[1]+er[1]), 'beta_lin_slope'=list(c(0.0)), 'beta_lin_int'=list(c(0.0)),'beta_quad'=list(c(0.0)),'beta_step_slope'=list(c(0.0)),'beta_step_int'=list(c(0.0)))
            e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir-1, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,length(uniq))
            xu <- e$x
            yu <- e$y
    #        stop()
            print(length(x))
            print("plotting risk")
            jpeg(paste("risk_plot_",fir,".jpg",sep=""), units="in", width=5, height=5, res=1200)
            plot(x,y, type="l", xlab=all_names[fir],ylab="Relative Risk",col='black',ylim=c(min(yl),max(yu)))
            lines(xl,yl,col='black')
            lines(xu,yu,col='black')
            dev.off()
            #
        } else {
            #
            lfir <- c(all_names[fir])
            uniq <- unlist(unique(df[,..lfir]))
    #        print(uniq)
            #
            e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir-1, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,length(uniq))
            x <- e$x
            y <- e$y
            #
            a_loglin[fir-1] = b[fir]-er[fir]
            e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir-1, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,length(uniq))
            xl <- e$x
            yl <- e$y
            #
            a_loglin[fir-1] = b[fir]+er[fir]
            e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir-1, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu,Plot_Type,length(uniq))
            xu <- e$x
            yu <- e$y
    #        stop()
            print(length(x))
            print("plotting risk")
            jpeg(paste("risk_plot_",fir,".jpg",sep=""), units="in", width=5, height=5, res=1200)
            plot(x,y, type="p", xlab=all_names[fir],ylab="Relative Risk",col='black',ylim=c(min(yl),max(yu)))
            lines(xl,yl,col='black',type="b")
            lines(xu,yu,col='black',type="b")
            dev.off()
            #
        }
        #
    }
}

time_var_coxph <-function(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event){
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
    print(length(tu))
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
    ce <- c(time1,time2,event)
    x_lin=as.matrix(df[,..lin_n])
    x_loglin=as.matrix(df[,..loglin_n])
    x_plin=as.matrix(df[,..plin_n])
    x_dose=as.matrix(df[,..dose_n])
    print("all results")
    e <- peanut_schoenfeld_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df[,..ce]),tu)
    for (i in 1:ncol(e)){
        print(i)
        pname<-paste("rplot",i,".jpg",sep="")
        df_res <- data.frame(x=tu,y=e[,i])
        jpeg(pname, units="in", width=5, height=5, res=1200)
        smoothScatter(df_res$x,df_res$y, xlab="age",ylab="Residual",main=all_names[i],nbin=200, colramp= colorRampPalette(c("white","red")))
        #
        level=0.95
        smodel <- loess(y ~ x, data=df_res)
        pred <- predict(smodel, se = TRUE)
        y = pred$fit
        ci <- pred$se.fit * qt(level / 2 + .5, pred$df)
           #
        lines(y - ci ~ tu,lwd=2.0)
        lines(y + ci ~ tu,lwd=2.0)
        lines(y ~ tu,lwd=2.0)
        dev.off()
        #
    }
}

Stratified_Baseline <- function(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event,strat_cov,strat_count){
    #-------------------------------------------------------------------------------------------------------------#
    #   df is the data that will be used
    #   df is changed to a data.table to make filtering easier
    colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)
    df <- df[sexm==1,]
    setkeyv(df, c(time2, event,time1))
    #-------------------------------------------------------------------------------------------------------------#
    #   The goal is to precompute the indices needed for each event time
    #   The file has structure {start of all data, end of all data, first event index, second event index, ...,}
    #
    #
    #
    if (length(lin_n)!=0){
        lin_n = Check_Dupe_Columns(df,lin_n)
    }
#    print(lin_n)
    if (length(loglin_n)!=0){
        loglin_n = Check_Dupe_Columns(df,loglin_n)
    }
    if (length(plin_n)!=0){
        plin_n = Check_Dupe_Columns(df,plin_n)
    }
    if (length(dose_paras$names)!=0){
        dose_paras$names = Check_Dupe_Columns(df,dose_paras$names)
    }
    #
    dose_n <- dose_paras$names
    #
    all_names <- c()
    all_names <- c(all_names, dose_paras$terms)
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
    ce <- c(time1,time2,event)
    #
    dfend <- df[get(event)==1, ]
    cov_vals <- unlist(unique(dfend[,..strat_cov]))
    if (length(cov_vals)>strat_count){
        strat_count <- length(cov_vals)
    }
    #
    Prob_vals <- quantile(cov_vals,probs=seq(0,1,1/strat_count),names=FALSE)
    print(Prob_vals)
    for (check_i in 1:strat_count){
        df_0 <- df[(get(strat_cov)>=Prob_vals[check_i])&(get(strat_cov)<=Prob_vals[check_i+1]),]
        x_lin=as.matrix(df_0[,..lin_n])
        x_loglin=as.matrix(df_0[,..loglin_n])
        x_plin=as.matrix(df_0[,..plin_n])
        x_dose=as.matrix(df_0[,..dose_n])
        dfend <- df_0[get(event)==1, ]
        #
        tu <- unlist(unique(dfend[,..time2]))
        print(length(tu))
        ##
    #    print(length(tu))
        ce <- c(time1,time2,event)
        e <- peanut_transition(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df_0[,..ce]),tu)
        #
        a_lin <- e$Parameter_Lists$beta_lin
        a_loglin <- e$Parameter_Lists$beta_loglin
        a_plin <- e$Parameter_Lists$beta_plin
        dose_paras$beta_loglin_top <- e$Parameter_Lists$beta_loglin_tops
        #
        Plot_Type=c("SURV","not_used")
        #                      a_lin,        a_loglin,        a_plin,   x_lin,  x_loglin,  x_plin, x_dose, fir, modelform, ntime,    include_bool,  Control,  Dose_paras,  d  f_groups,  tu 
        e <- peanut_plot(c(0.0,a_lin), c(0.0,a_loglin), c(0.0,a_plin),  x_lin,  x_loglin,  x_plin, x_dose, fir, modelform,length(tu),term_bool, control, dose_paras,as.matrix(df_0[,..ce]),tu,Plot_Type,0)
        t <- c()
    #    ch <- c()
        surv <- c()
        dt <- 1
        dft=data.table("time"=tu,"base"=e$baseline)
        for (i in tu[1]:tu[length(tu)]){
            t <- c(t,i)
            temp <- sum(dft[time<i, base])
    #        ch <- c(ch, temp)
            surv <- c(surv, exp(-1*temp))
        }
        jpeg(paste('surv_',strat_cov,'_plot_',check_i,'.jpg',sep=""), units="in", width=5, height=5, res=1200)
        plot(t,surv, type="s", xlab="age",ylab="Survival")
        dev.off()
    }
}

#fname <- '../Combined MW NPP IR Lung Lag10 - TEST 4.18.22.csv'
fname <- '/home/user/Documents/CMS/Combined_MW_NPP_IR_Lung_Lag10_-_TEST_4.18.22.csv'
lin_n <- c()
loglin_n <- c("sexm","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
plin_n <- c()

a_lin=c()
a_loglin <- c(-0.0006300516,-0.0630181076,-0.3932366917,-0.6613753288,-1.1079753794,-0.7946066410,-0.1990600045,0.1616788011,-0.6317184838,-0.0572978471,-0.5560610214,-1.1076729774,-1.2444050074,-1.4178535299)
a_plin=c()
a_dose=c(0.0001946145)

modelform <- 'M'
fir=0

control <- list('lr' = 0.75,'maxiter' = 1,'halfmax' = 1,'epsilon' = 1e-5,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0)

time1="age_entry"
time2="age_exit"
event="lung"

dose_paras <- list('names' =c('cumulative_dose_lung_lag10'),'terms'=c('beta_loglin_top'), 'beta_loglin_slope'=list(c(1.0)), 'beta_loglin_top'=list(c(0.0001946145)), 'beta_lin_slope'=list(c(0.0)), 'beta_lin_int'=list(c(0.0)),'beta_quad'=list(c(0.0)),'beta_step_slope'=list(c(0.0)),'beta_step_int'=list(c(0.0)))
#
plot_coxph(fname,lin_n,loglin_n,plin_n,a_lin,a_loglin,a_plin,modelform,dose_paras,fir,control,time1,time2,event)


#
lin_n <- c()
loglin_n <- c("sexm","YOB1","YOB2","YOB3","YOB4","COH_EDUC1","COH_EDUC2","COH_EDUC3","COH_EDUC4","COH_EDUC5","COH_EDUC6","COH_EDUC7","COH_EDUC8","COH_EDUC9")
plin_n <- c()

a_lin=c()
a_loglin <- rep(-.01,length(loglin_n)-1)
a_plin=c()
a_dose=c(-0.1)

modelform <- 'M'
fir=0

control <- list('lr' = 0.75,'maxiter' = 30,'halfmax' = 5,'epsilon' = 1e-7,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-7, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=10.0)

for (i in 1:length(a_loglin)){
    cov_name <- loglin_n[i]
    Stratified_Baseline(fname,lin_n,loglin_n,plin_n,a_lin,setdiff(a_loglin,c(cov_name)),a_plin,modelform,dose_paras,fir,control,time1,time2,event,cov_name,2)
}




