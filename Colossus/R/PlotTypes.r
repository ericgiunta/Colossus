#' calculates and plots martingale residuals with a named dose column
#' \code{CoxMartingale} uses user provided data, columns, and identifier to create plots
#'
#' @param verbose boolean controlling additional printing
#' @param df data with covariates and times
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param e output from a baseline calculation
#' @param t event times
#' @param ch cumulative hazards
#' @param dnames list of covariate columns to plot by
#' @param Plot_Name plot identifier
#' @param age_unit age unit
#' @param studyID id to group by, NaN for no grouping
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
#'
CoxMartingale <- function(verbose, df, time1, time2, event,e, t, ch, dnames, Plot_Name, age_unit, studyID){
    IDS <- base <- res <- doses <- NULL
    if (verbose){
        print("Martingale Plots")
    }
    #
    time_s <- df[,get(time1)]
    time_e <- df[,get(time2)]
    ch_fun <- approxfun(x=t,y=ch,rule=2)
    ch_e <- ch_fun(time_e)
    ch_s <- ch_fun(time_s)
    #
    e_i <- df[,get(event)]
    if (verbose){
        print("Cumulative Hazard Functions Estimated")
    }
    for (cov_i in 1:length(dnames)){
        dname <- dnames[cov_i]
        if (verbose){
            print(paste("Martingale Plot: ",dname,sep=""))
        }
        if (studyID %in% names(df)){
            dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"IDS"=unlist(df[,studyID,with=FALSE],use.names=FALSE),"time"=time_e,"cov"=unlist(df[, dname, with = FALSE],use.names=FALSE))
            #
            name_temp <- names(dfr)
            for (i in 1:length(name_temp)){
                if (grepl( "cov", name_temp[i], fixed = TRUE)){
                    setnames(dfr,name_temp[i],c("cov"))
                } else if (grepl( "IDS", name_temp[i], fixed = TRUE)){
                    setnames(dfr,name_temp[i],c("IDS"))
                }
            }
            #
            dfr$res = dfr$e - (dfr$ch_e-dfr$ch_s) * dfr$Risks
            #
            #
            Martingale_Error <- dfr[, lapply(.SD,sum), by=IDS]
            times <- dfr[, lapply(.SD,max), by=IDS]
        } else {
            dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"time"=time_e,"cov"=unlist(df[, dname, with = FALSE],use.names=FALSE))
            #
            name_temp <- names(dfr)
            for (i in 1:length(name_temp)){
                if (grepl( "cov", name_temp[i], fixed = TRUE)){
                    setnames(dfr,name_temp[i],c("cov"))
                } else if (grepl( "IDS", name_temp[i], fixed = TRUE)){
                    setnames(dfr,name_temp[i],c("IDS"))
                }
            }
            #
            dfr$res = dfr$e - (dfr$ch_e) * dfr$Risks
            #
            Martingale_Error <- dfr
            times <- dfr
        }
        #
        ##
        if (verbose){
            print(paste("Martingale Plotting",sep=""))
        }
        dft <- data.table("cov_max"=times$cov,"time_max"=times$time,"res_sum"=Martingale_Error$res)
        ##
        g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$cov_max, y=.data$res_sum)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("Max",dname,sep=" "), y="Martingale Residuals")
        ggplot2::ggsave(paste('martin_plot_',dname,Plot_Name,'.jpg',sep="_"),device="jpeg",dpi="retina")
        ##
    }
    #
    if (studyID %in% names(df)){
        dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"IDS"=unlist(df[,studyID,with=FALSE],use.names=FALSE),"time"=time_e)
        #
        #
        name_temp <- names(dfr)
        for (i in 1:length(name_temp)){
            if (grepl( "IDS", name_temp[i], fixed = TRUE)){
                setnames(dfr,name_temp[i],c("IDS"))
            }
        }
        #
        dfr$res = dfr$e - (dfr$ch_e-dfr$ch_s) * dfr$Risks
        #
        Martingale_Error <- dfr[, lapply(.SD,sum), by=IDS]
        times <- dfr[, lapply(.SD,max), by=IDS]
    } else {
        dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"time"=time_e)
        #
        name_temp <- names(dfr)
        for (i in 1:length(name_temp)){
            if (grepl( "IDS", name_temp[i], fixed = TRUE)){
                setnames(dfr,name_temp[i],c("IDS"))
            }
        }
        #
        dfr$res = dfr$e - (dfr$ch_e) * dfr$Risks
        #
        Martingale_Error <- dfr
        times <- dfr
    }
    dft <- data.table("time_max"=times$time,"res_sum"=Martingale_Error$res)
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$time_max, y=.data$res_sum)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("Max Age",sep=" "), y="Martingale Residuals")
    ggplot2::ggsave(paste('martin_plot',Plot_Name,'.jpg',sep='_'),device="jpeg",dpi="retina")
    ##
    return ("passed")
}

#' calculates and plots survival plots of the estimated baseline
#' \code{CoxSurvival} uses user provided data, columns, and identifier to create plots
#'
#' @param t event times
#' @param h hazards of baseline
#' @param ch cumulative hazards of baseline
#' @param surv survival fraction of baseline
#' @param Plot_Name plot identifier
#' @param verbose boolean controlling additional printing
#' @param time_lims limits for x axis of plot
#' @param age_unit age unit
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
CoxSurvival <- function(t,h,ch,surv,Plot_Name,verbose,time_lims, age_unit){
    if (verbose){
        print("Survival Plots")
    }
    #
    dft <- data.table("t"=t,"h"=h,"ch"=ch,"surv"=surv)
    setkeyv(dft,"t")
    dft <- dft[(t>=time_lims[1])&(t<=time_lims[2]),]
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$t, y=.data$ch)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Cumulative Hazard")
    ggplot2::ggsave(paste('ggplot2_ch_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$t, y=.data$surv)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Survival")
    ggplot2::ggsave(paste('ggplot2_surv_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$t, y=.data$h)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Hazard Estimate")
    ggplot2::ggsave(paste('ggplot2_H_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    #
    return ("passed")
}

#' calculates and plots Kaplan-Mieier survival plots
#' \code{CoxKaplanMeier} uses user provided data, columns, and identifier to create plots, plots the survival for subjects above/below the mean column value for each column
#'
#' @param verbose boolean identifying if extra plotting information should be written to the console
#' @param verbosec boolean identifying if extra regression information should be written to the console           
#' @param studyID column identifying each individual                 
#' @param names columns names for elements of the model used to identify data columns      
#' @param df data used for regression              
#' @param event column used for event status             
#' @param time1 column used for time period starts            
#' @param time2 column used for time period end            
#' @param tu unique event times                 
#' @param Term_n term numbers for each element of the model          
#' @param tform subterm type for each element of the model          
#' @param a_n starting parameters for regression      
#' @param er standard deviation for the parameters        
#' @param fir term number for the initial term used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at only used for testing runs with a single varying parameter
#' @param modelform string specifying the model type             
#' @param control list of parameters controlling the convergence            
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant       
#' @param Plot_Type list of parameters controlling the plot options  
#' @param age_unit age unit
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
CoxKaplanMeier <- function(verbose, verbosec, studyID,names,df,event,time1,time2,tu,Term_n, tform, a_n, er, fir, der_iden, modelform, control,keep_constant, Plot_Type, age_unit){
    if (verbose){
        print("KM Plots")
    }
    base  <- NULL
    ce <- c(time1,time2,event)
    all_names <- unique(names)
    dfc <- match(names,all_names)
    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
#    stop()
    #
    df_study <- df[, lapply(.SD,max), by=studyID]
    #
    for (fir_KM in 1:length(all_names)){
        if (verbose){
            print(fir_KM)
        }
        dfend <- df_study[get(event)==1, ]
        fir_c <- mean(dfend[,get(all_names[fir_KM])])
        fmean <- mean(fir_c) #average value to split by
        df_u <- df_study[get(all_names[fir_KM])>=fmean,] #data above the mean
        df_l <- df_study[get(all_names[fir_KM])<=fmean,] #data below the mean
        #
        df_u_end <- df_u[get(event)==1, ]
        df_l_end <- df_l[get(event)==1, ]
        #
        u_num = length(unlist(unique(df_u[,studyID, with = FALSE]),use.names=FALSE)) #number of unique individuals in upper set
        l_num = length(unlist(unique(df_l[,studyID, with = FALSE]),use.names=FALSE)) #number of unique individuals in lower set
        t_num = length(unlist(unique(df_study[,studyID, with = FALSE]),use.names=FALSE)) #number of unique individuals in total set
        #
        t_u <- c(0) #time data
        t_l <- c(0)
        t_t <- c(0)
        n_u <- c(1) #surviving decimal
        n_l <- c(1)
        n_t <- c(1)
        iden <- c("above","below","combined") #list of which set
#        n_l <- c(1)
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)) #all event times
        #
        if (verbose){
            print("Starting KM Plots")
        }
#        stop()
        #
        tu_s <- c(0.0,tu)
        for (i_t in 1:length(tu)){
            i <- tu[i_t]
            #
#            print("----------")
#            print(i)
            #
            df0 <- df_u_end[get(time2)<=i,] #set of all intervals prior to this point in upper set
            df0 <- df0[(get(time2)>tu_s[i_t]),]
            df1 <- df_u[(get(time2)>tu_s[i_t]),]
            u_ev = sum(df0[, get(event)]) #number of intervals with event in upper set prior to the time point
            u_num = nrow(df1)
#            print("0")
            df0 <- df_l_end[get(time2)<=i,] #set of all intervals prior to this point in lower set
            df0 <- df0[(get(time2)>tu_s[i_t]),]
            df1 <- df_l[(get(time2)>tu_s[i_t]),]
            l_ev = sum(df0[, get(event)]) #number of intervals with event in lower set prior to the time point
            l_num = nrow(df1)
#            print("1")
            df0 <- dfend[get(time2)<=i,] #set of all intervals prior to this point in lower setv
            df0 <- df0[(get(time2)>tu_s[i_t]),]
            df1 <- df_study[(get(time2)>tu_s[i_t]),]
            t_ev = sum(df0[, get(event)]) #number of intervals with event in lower set prior to the time point
            t_num = nrow(df1)
#            print("2")
            #
            if (u_ev>0){ #if there are events in the upper set
                t_u <- c(t_u,i) #adds new time
                temp <- (u_num - u_ev)/u_num #compares number of events to maximum total
                n_u <- c(n_u, n_u[length(n_u)]*temp) #adds to list of survival decimals
            }
            if (l_ev>0){
                temp <- (l_num - l_ev)/l_num
                t_l <- c(t_l,i)
                n_l <- c(n_l, n_l[length(n_l)]*temp)
            }
            if (t_ev>0){
                temp <- (t_num - t_ev)/t_num
                t_t <- c(t_t,i)
                n_t <- c(n_t, n_t[length(n_t)]*temp)
            }
#            print("3")
        }
#        print(iden)
        n <- c(n_u,n_l,n_t)
        t <- c(t_u,t_l,t_t)
        iden <- c(rep("above",length(n_u)),rep("below",length(n_l)),rep("combined",length(n_t)))
        #
        dft <- data.table("t_u"=t,"n_u"=n,"iden"=iden)
        #
        g <- ggplot2::ggplot(data=dft,ggplot2::aes(x=.data$t_u, y=.data$n_u,color=.data$iden)) + ggplot2::geom_line() + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Survival")
        g <- g + ggplot2::scale_color_discrete(name = all_names[fir_KM])
        ggplot2::ggsave(paste("KM_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
#        stop()
        #
        df_u <- df[get(all_names[fir_KM])>=fmean,] #data above the mean
        df_l <- df[get(all_names[fir_KM])<=fmean,] #data below the mean
        #
        x_all=as.matrix(df_l[,all_names, with = FALSE])
        dfend <- df_l[get(event)==1, ]
        #
        #
        if (verbose){
            print("Lower Survival Plots")
        }
        #
        #
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        e0 <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_l[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        tu0 <- tu
        ##
        x_all=as.matrix(df_u[,all_names, with = FALSE])
        dfend <- df_u[get(event)==1, ]
        #
        #
        if (verbose){
            print("Upper Survival")
        }
        #
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        e1 <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_u[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        tu1 <- tu
        ##
        x_all=as.matrix(df[,all_names, with = FALSE])
        dfend <- df[get(event)==1, ]
        #
        #
        if (verbose){
            print("Total Survival")
        }
        #
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        e2 <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        tu2 <- tu
        tmin <- min(c(min(tu0),min(tu1),min(tu2)))
        t <- c(tmin,tmin,tmin)
        surv <- c(1,1,1)
        iden <- c('below','above','combined')
        dt <- 1
        if (verbose){
            print("writing survival")
        }
        dft0=data.table("time"=tu0,"base"=e0$baseline)
        dft1=data.table("time"=tu1,"base"=e1$baseline)
        dft2=data.table("time"=tu2,"base"=e2$baseline)
        dfend <- df[get(event)==1, ]
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
        for (i in tu[1]:tu[length(tu)]){
            t <- c(t,i)
            t <- c(t,i)
            t <- c(t,i)
            temp0 <- sum(dft0[time<i, base])
            temp1 <- sum(dft1[time<i, base])
            temp2 <- sum(dft2[time<i, base])
            surv <- c(surv, exp(-1*temp0))
            surv <- c(surv, exp(-1*temp1))
            surv <- c(surv, exp(-1*temp2))
            iden <- c(iden,"below","above","combined")
        }
        #
        Ls <- log(surv)
        Lls_u <- log(-Ls)
        Lt_u <- log(t)
        ##
        #
        if (verbose){
            print("Split log-log Survival Plots")
        }
        #
        #
        dft <- data.table("t"=Lt_u,"s"=Lls_u,"iden"=iden)
        #

        g <- ggplot2::ggplot(data=dft,ggplot2::aes(x=.data$t, y=.data$s,color=.data$iden)) + ggplot2::geom_line() + ggplot2::labs(x="Log-Age", y="Log of Log Survival")
        g <- g + ggplot2::scale_color_discrete(name = all_names[fir_KM]) 
        ggplot2::ggsave(paste("log_log_surv_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
    }
    ;
    return ("passed")
}


#' calculates and plots relative risk
#' \code{CoxRisk} uses user provided data, columns, and identifier to create plots of risk by covariate value for each column
#'
#' @param verbose boolean to control if additional information is printed to the console
#' @param df data used for regression              
#' @param event column used for event status             
#' @param time1 column used for time period starts            
#' @param time2 column used for time period end                  
#' @param names columns names for elements of the model used to identify data columns               
#' @param Term_n term numbers for each element of the model          
#' @param tform subterm type for each element of the model          
#' @param a_n starting parameters for regression              
#' @param fir term number for the initial term used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at only used for testing runs with a single varying parameter
#' @param modelform string specifying the model type             
#' @param control list of parameters controlling the convergence            
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant       
#' @param Plot_Type list of parameters controlling the plot options
#' @param b optimum parameter values used
#' @param er standard deviation of optimum parameters used
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
CoxRisk <- function(verbose,df, event, time1, time2, names,Term_n, tform, a_n, fir, der_iden, modelform, control,keep_constant, Plot_Type, b, er){
    fir_KM=0
    dfend <- df[get(event)==1, ]
    #
    ce <- c(time1,time2,event)
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    #
    tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
    if (verbose){
        print("start risk calculations")
    }
    for (fir_KM in 1:length(names)){
        lfir <- c(names[fir_KM])
        uniq <- unlist(unique(df[,lfir, with = FALSE]), use.names=FALSE)
        #
        der_iden <- fir_KM-1
#        print(lfir)
#        print(min(uniq))
#        print(max(uniq))
        if (TRUE){
            #    cox_ph_plot(Term_n, tform, a_n, a_er,dfc,x_all, fir, der_iden, modelform, Control, df_groups,                        tu, KeepConstant,  term_tot, Plot_Type , uniq_v)
            #
#            print("in")
            e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type ,length(uniq))
#            print("out")
            x <- e$x
            y <- e$y
            #
#            stop()
#            print(length(x))
#            print(length(y))
            if (TRUE){
                dft <- data.table("x"=x,"y"=y)
#                print(c(x[1],x[length(x)]))
                if (length(uniq)>100){
                    #
                    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$x, y=.data$y)) + ggplot2::geom_line(color="black") + ggplot2::labs(x=names[fir_KM], y="Relative Risk")
                    ggplot2::ggsave(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
                    #
                } else {
                    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$x, y=.data$y)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=names[fir_KM], y="Relative Risk")
                    ggplot2::ggsave(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
                    #
                }
            } else {
#                print(length(x))
                a_n[fir_KM] <- b[fir_KM] - er[fir_KM]
                #
                e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , length(uniq))
                xl <- e$x
                yl <- e$y
                #
                a_n[fir_KM] <- b[fir_KM] + er[fir_KM]
                #
                e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , length(uniq))
                xu <- e$x
                yu <- e$y
                #
                dft <- data.table("xu"=xu,"yu"=yu,"xl"=xl,"yl"=yl,"x"=x,"y"=y)
                if (length(uniq)>100){
                    #
                    print(length(x))
                    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$x, y=.data$y)) + ggplot2::geom_line(color="black") + ggplot2::labs(x=names[fir_KM], y="Relative Risk")
                    g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xl, y=.data$yl), color="black")
                    g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xu, y=.data$yu), color="black")
                    ggplot2::ggsave(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
                    #
                } else {
                    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$x, y=.data$y)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=names[fir_KM], y="Relative Risk")
                    g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xl, y=.data$yl), color="black")
                    g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xu, y=.data$yu), color="black")
                    ggplot2::ggsave(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
                    #
                }
            }
        }
        #
    }
    ;
    return ("passed")
}


#' calculates and plots survival curves for each unique value of the stratification column
#' \code{CoxStratifiedSurvival} uses user provided data, columns, and identifier to calculate the survival fraction for each strata
#'
#' @param verbose boolean to control if additional information is printed to the console
#' @param df data used for regression              
#' @param event column used for event status             
#' @param time1 column used for time period starts            
#' @param time2 column used for time period end                  
#' @param names columns names for elements of the model used to identify data columns               
#' @param Term_n term numbers for each element of the model          
#' @param tform subterm type for each element of the model          
#' @param a_n starting parameters for regression       
#' @param er standard deviation for the parameters               
#' @param fir term number for the initial term used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at only used for testing runs with a single varying parameter
#' @param modelform string specifying the model type             
#' @param control list of parameters controlling the convergence            
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant       
#' @param Plot_Type list of parameters controlling the plot options
#' @param Strat_Col column to stratify by
#' @param time_lims limits for x axis of plot
#' @param age_unit age unit
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
CoxStratifiedSurvival <- function(verbose, df, event, time1, time2, names,Term_n, tform, a_n, er, fir, der_iden, modelform, control,keep_constant, Plot_Type, Strat_Col,time_lims, age_unit){
    setkeyv(df, c(time2, event, Strat_Col))
    dfend <- df[get(event)==1, ]
    base  <- NULL
    #
    ce <- c(time1,time2,event,Strat_Col)
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    #
    tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
    uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
    e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
    a_n <- e$beta_0
    Plot_Name <- Plot_Type[2]
    if (verbose){
        print(all_names)
        print(uniq)
        print(a_n)
    }
#    jpeg(paste('surv_plot',Strat_Col,Plot_Name,'.jpg',sep="_"), units="in", width=5, height=5, res=1200)
    tt <- c()
    tsurv <- c()
    categ <- c()
    for (col_u in uniq){
        if (verbose){
            print(paste("Starting Stratification: ",col_u))
        }
        df0 <- df[get(Strat_Col)==col_u,]
        dfend <- df0[get(event)==1, ]
        x_all=as.matrix(df0[,all_names, with = FALSE])
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        #
        e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df0[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        #
        t <- c()
        h <- c()
        ch <- c()
        surv <- c()
        dt <- 1
        if (verbose){
            print("writing survival data")
#            print(time_lims)
        }
        dft=data.table("time"=tu,"base"=e$baseline)
        if (verbose){
            print(dft)
            ;
        }
        for (i in tu){
            if ((i<=time_lims[2])&(i>=time_lims[1])){
                t <- c(t,i)
                temp <- sum(dft[time<i, base])
                ch <- c(ch, temp)
                if (length(h)==0){
                    h <- c(temp)
                } else {
                    h <- c(h, ch[length(ch)]-ch[length(ch)])
                }
                surv <- c(surv, exp(-1*temp))
            }
        }
        tt <- c(tt, t)
        tsurv <- c(tsurv, surv)
        categ <- c(categ, rep(paste(col_u),length(t)))
    }
    if (verbose){
        print("combining data")
    }
    dft=data.table("t"=tt,"surv"=tsurv,"cat_group"=categ)
    sbreaks <- c()
    slabels <- c()
    for (i in 1:length(uniq)){
        sbreaks <- c(sbreaks, paste(uniq[i]))
        slabels <- c(slabels, paste("For ",Strat_Col,"=",uniq[i],sep=""))
    }
    if (verbose){
        print("plotting survival data")
#            print(t)
    }
    g <- ggplot2::ggplot() + ggplot2::geom_point(data=dft,ggplot2::aes(x=.data$t, y=.data$surv,group=.data$cat_group,color=.data$cat_group))#, direction="vh")
    g <- g + ggplot2::scale_colour_discrete(breaks=sbreaks, labels=slabels)
    g <- g + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Survival") + ggplot2::ylim(0,1)
    if (verbose){
        print("saving survival data")
    }
    ggplot2::ggsave(paste('strat_surv_plot',Strat_Col,Plot_Name,'.jpg',sep="_"),device="jpeg",dpi="retina")
    ;
    return ("passed")
}
        
        
#' calculates a smoothing factor
#' \code{K_Smooth} calculates 3/4 of (1-x)^2
#'
#' @param x point to smooth using
#'
#' @return returns the smoothing factor
#' @export
#'
K_Smooth <- function(x){
    return (0.75*(1-x^2))
}

#' calculates and plots survival plots of the smoothed estimated baseline hazard
#' \code{CoxSmoothHazard} uses user provided data, band-width, and identifier to create plots
#'
#' @param dft dataframe of times and hazards
#' @param Plot_Name plot identifier
#' @param verbose boolean controlling additional printing
#' @param bw band-width used for smoothing
#' @param time_lims limits for x axis of plot
#' @param age_unit age unit
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
CoxSmoothHazard <- function(dft,Plot_Name,verbose,bw,time_lims, age_unit){ 
    if (verbose){
        print("writing survival data")
    }
    tu=dft$time
    base=dft$base
    base_er=dft$basehaz
    ti <- seq(time_lims[1], time_lims[2], length.out = 1000)
    h <- rep(0.0,length(ti))
    h_l <- rep(0.0,length(ti))
    h_u <- rep(0.0,length(ti))
    er_h <- rep(0.0,length(ti))
    for (j in 1:length(ti)){
        for (i in 1:length(tu)){
            if (abs(ti[j]-tu[i])<=bw){
                xtemp <- (ti[j]-tu[i])/bw
                h[j]   <- h[j]   + (base[i]*K_Smooth(xtemp))
                er_h[j] <- er_h[j] + (base_er[i] * K_Smooth(xtemp)^2)
#                h_l[j] <- h_l[j] + (base[i]-sqrt(base_er[i]))*K_Smooth(xtemp)
#                h_u[j] <- h_u[j] + (base[i]+sqrt(base_er[i]))*K_Smooth(xtemp)
            }
        }
    }
    for (j in 1:length(ti)){
        h[j] = h[j] / bw
        er_h[j] = sqrt(er_h[j]) / bw
#        print(er_h[j])
        h_u[j] = h[j] + er_h[j]
        h_l[j] = h[j] - er_h[j]
    }
    dft <- data.table("ti"=ti,"h"=h,"h_l"=h_l,"h_u"=h_u)
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$ti, y=.data$h)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Survival")
#    g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$ti, y=.data$h_u), color="black")
#    g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$ti, y=.data$h_l), color="black")
    ggplot2::ggsave(paste('H_smooth_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    ;
    return ("passed")
}
        
        
#' Calculates Schoenfeld residuals for a Cox Proportional Hazards regression and plots
#' \code{RunCox_Schoenfeld_Residual} uses user provided data, time/event columns, vectors specifying the model, and options to calculate the residuals
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param age_unit age unit
#' @param Plot_Name plot identifier
#' @param alpha significance level for two tail t test
#'
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data

PlotCox_Schoenfeld_Residual <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,age_unit,Plot_Name,alpha=0.05){        
    setkeyv(df, c(time2, event))
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(a_n)-length(names)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    #
    a_n0 <- copy(a_n)
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    res_list <- cox_ph_schoenfeld_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    res <- res_list$residual
    res_scaled <- res_list$scaled
    degree_f <- res_list$df
    #
    for (cov in 1:length(a_n)){
        print(names[cov])
        y <- unlist(res[,cov],use.names=FALSE)
        y_scale <- unlist(res_scaled[,cov],use.names=FALSE)
        #
        dft <- data.table("time"=tu,"y"=y)
        dft$y_scale <- y_scale
        #
        g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$time, y=.data$y)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y=paste("Schoenfeld Residual (",names[cov],")",sep=" "))
        ggplot2::ggsave(paste("schoenfeld_",cov,"_",Plot_Name,".jpg",sep=""),device="jpeg",dpi="retina")
        #
        g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$time, y=.data$y_scale)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y=paste("Schoenfeld Residual Scaled (",names[cov],")",sep=" "))
        ggplot2::ggsave(paste("schoenfeld_scaled_",cov,"_",Plot_Name,".jpg",sep=""),device="jpeg",dpi="retina")
        #
    }
    return ("Passed")
}
          
        
#' Calculates and returns data for time by cumulative hazard
#' \code{GetCensWeight} uses user provided data, time/event columns, vectors specifying the model, and options generate an estimate of the censoring rate, plots, and returns the data
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#' @param plot_options list of parameters controlling the plot options
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
#'
GetCensWeight <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options){
    if (plot_options$verbose){
        print("Starting Plot Function")
    }
    setkeyv(df, c(time2, event))
    base  <- NULL
    der_iden <- 0
    Plot_Name <- plot_options$name
    Plot_Type <- "SURV"
    if (plot_options$verbose){
        print("Getting Plot Info")
    }
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (plot_options$verbose){
            print("no events")
        }
        stop()
    }
    if (plot_options$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    #
    if ("age_unit" %in% names(plot_options)){
        ;
    } else {
        plot_options$age_unit <- "unitless"
    }
    time_lims <- c(min(tu),max(tu))
    for (iden_col in c("verbose")){
        if (iden_col %in% names(plot_options)){
            ;
        } else {
            plot_options[iden_col] <- FALSE
        }
    }
    #
    control <- Def_Control(control)
    verbose <- copy(plot_options$verbose)
    verbosec <- copy(control$verbose)
    maxiterc <- copy(control$maxiter)
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    all_names <- unique(names)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    time1 <- ce[1]
    time2 <- ce[2]
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    control$maxiter <- -1
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    control$maxiter <- maxiterc
    b <- e$beta_0
    er <- e$Standard_Deviation
    #
    if (verbose){
        print("starting ph_plot")
    }
    #
    e <- cox_ph_plot(Term_n, tform, a_n,er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
    #
    t <- c()
    h <- c()
    ch <- c()
    surv <- c()
    dt <- 1
    if (verbose){
        print("writing survival data")
    }
    dft=data.table("time"=tu,"base"=e$baseline,"basehaz"=e$standard_error)
    for (i in tu){
        t <- c(t,i)
        temp <- sum(dft[time<i, base])
        ch <- c(ch, temp)
        if (length(h)==0){
            h <- c(temp)
        } else {
            h <- c(h, ch[length(ch)]-ch[length(ch)-1])
        }
        surv <- c(surv, exp(-1*temp))
    }
    if (verbose){
        print("Survival Plots")
    }
    age_unit <- plot_options$age_unit
    #
    dft <- data.table("t"=t,"h"=h,"ch"=ch,"surv"=surv)
    setkeyv(dft,"t")
    dft <- dft[(t>=time_lims[1])&(t<=time_lims[2]),]
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$t, y=.data$ch)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Cumulative Hazard")
    ggplot2::ggsave(paste('weight_ch_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$t, y=.data$surv)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Survival")
    ggplot2::ggsave(paste('weight_surv_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    #
    g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$t, y=.data$h)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Hazard Estimate")
    ggplot2::ggsave(paste('weight_H_plot',Plot_Name,".jpeg",sep="_"),device="jpeg",dpi="retina")
    return (dft)
}
        
        
        
