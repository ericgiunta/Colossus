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
#' @param dname dose column
#' @param Plot_Name plot identifier
#' @param age_unit age unit
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
#'
CoxMartingale <- function(verbose, df, time1, time2, event,e, t, ch, dname, Plot_Name, age_unit){
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
    dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"IDS"=df$studyID,"time"=time_e,"doses"=df[, dname, with = FALSE])
    dfr$res = dfr$e - (dfr$ch_e-dfr$ch_s) * dfr$Risks
    #
    Martingale_Error <- dfr[, lapply(.SD,sum), by=IDS]
    times <- dfr[, lapply(.SD,max), by=IDS]
    ##
    jpeg(paste('dose_sum_by_max',Plot_Name,'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
    plot(times[,doses],Martingale_Error[,doses], xlab="Max Dose",ylab="Sum of Doses",type='p')
    dev.off()
    ##
    jpeg(paste('martin_plot',Plot_Name,'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
    smoothScatter(times[,time],Martingale_Error[,res], xlab=paste("age (",age_unit,")",sep=""),ylab="Martingale Residuals",nbin=100, colramp= colorRampPalette(c("white","red")))
    dev.off()
    jpeg(paste('martin_plot_dose',Plot_Name,'.jpg',sep="_"), units="in", width=5, height=5, res=1200)
    smoothScatter(times[doses>0,doses],Martingale_Error[doses>0,res], xlab="Dose",ylab="Martingale Residuals",nbin=200, colramp= colorRampPalette(c("white","red")))
    dev.off()
    ;
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
#' @param fir term number for the intial term used for models of the form T0*f(Ti) in which the order matters
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
    #
    for (fir_KM in 1:length(all_names)){
        if (verbose){
            print(fir_KM)
        }
        dfend <- df[get(event)==1, ]
        fir_c <- mean(dfend[,get(all_names[fir_KM])])
        fmean <- mean(fir_c)
        df_u <- df[get(all_names[fir_KM])>=fmean,]
        df_l <- df[get(all_names[fir_KM])<=fmean,]
        #
        u_num = length(unlist(unique(df_u[,studyID, with = FALSE]),use.names=FALSE))
        l_num = length(unlist(unique(df_l[,studyID, with = FALSE]),use.names=FALSE))
        #
        t_u <- c(0)
        t_l <- c(0)
        n_u <- c(1)
        n_l <- c(1)
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        #
        if (verbose){
            print("Starting times Plots")
        }
        #
        for (i in tu[1]:tu[length(tu)]){
            #
            df0 <- df_u[(get(time1)<i)&(get(time2)<=i),]
            u_ev = sum(df0[, get(event)])
            df0 <- df_l[(get(time1)<i)&(get(time2)<=i),]
            l_ev = sum(df0[, get(event)])
            #
            if (u_ev>0){
                t_u <- c(t_u,i)
                temp <- (u_num - u_ev)/u_num
                n_u <- c(n_u, temp)
            }
            if (l_ev>0){
                t_l <- c(t_l,i)
                temp <- (l_num - l_ev)/l_num
                n_l <- c(n_l, temp)
            }
        }
#        jpeg(paste("KM_",fir_KM,"_",Plot_Type[2],".jpg",sep=""), units="in", width=5, height=5, res=1200)
#        plot(t_u,n_u, xlab=all_names[fir_KM],ylab="Survival Fraction",ylim=c(min(n_l),max(n_u)),col='red',type='l')
#        lines(t_l,n_l,col='blue')
#        legend("topright", legend=c("Below Mean", "Above Mean"), col=c("red", "blue"), lty=1:2, cex=0.8,title=all_names[fir_KM])
#        dev.off()
        #
        dft <- data.table("t_u"=t_u,"t_l"=t_l,"n_u"=n_u,"n_l"=n_l)
        #
        g <- ggplot2::ggplot() + ggplot2::geom_line(data=dft,ggplot2::aes(x=.data$t_u, y=.data$n_u, color="red"))+ ggplot2::geom_line(data=dft,ggplot2::aes(x=.data$t_u, y=.data$n_u), color="red") + ggplot2::labs(x=paste("age (",age_unit,")",sep=""), y="Survival")
        g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$t_l, y=.data$n_l, color="blue")) + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$t_l, y=.data$n_l), color="blue")
        g <- g + ggplot2::scale_color_manual(name = all_names[fir_KM], values = c("Above Mean" = "red", "Below Mean" = "blue"))      
        ggplot2::ggsave(paste("KM_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
        #
        if (verbose){
            print("Starting times Plots")
        }
        #
        x_all=as.matrix(df_l[,all_names, with = FALSE])
        dfend <- df_l[get(event)==1, ]
        #
        #
        if (verbose){
            print("Lower Survival Plots")
        }
        #
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_l[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        t <- c()
        ch <- c()
        surv <- c()
        dt <- 1
        if (verbose){
            print("writing survival")
        }
        dft=data.table("time"=tu,"base"=e$baseline)
        for (i in tu[1]:tu[length(tu)]){
            t <- c(t,i)
            temp <- sum(dft[time<=i, base])
            ch <- c(ch, temp)
            surv <- c(surv, exp(-1*temp))
        }
        #
        Ls <- log(surv)
        Lls_l <- log(-Ls)
        Lt_l <- log(t)
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
        e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_u[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        t <- c()
        ch <- c()
        surv <- c()
        dt <- 1
        if (verbose){
            print("writing survival")
        }
        dft=data.table("time"=tu,"base"=e$baseline)
        for (i in tu[1]:tu[length(tu)]){
            t <- c(t,i)
            temp <- sum(dft[time<i, base])
            ch <- c(ch, temp)
            surv <- c(surv, exp(-1*temp))
        }
        #
        Ls <- log(surv)
        Lls_u <- log(-Ls)
        Lt_u <- log(t)
        ##
        #
        if (verbose){
            print("Survival Plots")
        }
        #
#        jpeg(paste('log_log_surv_plot',fir_KM,"_",Plot_Type[2],'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
#        plot(Lt_l,LLs_l, xlab="Log-Age",ylab="Log of Log Survival",col='red',type='l')
#        lines(Lt_u,LLs_u,col='blue')
#        legend("topleft", legend=c("Below Mean", "Above Mean"), col=c("red", "blue"), lty=1:2, cex=0.8,title=all_names[fir_KM])
#        dev.off()
        #
        dft <- data.table("Lt_u"=Lt_u,"Lt_l"=Lt_l,"Lls_u"=Lls_u,"Lls_l"=Lls_l)
        #
        g <- ggplot2::ggplot() + ggplot2::geom_line(data=dft,ggplot2::aes(x=.data$Lt_u, y=.data$Lls_u, color="red")) + ggplot2::geom_line(data=dft,ggplot2::aes(x=.data$Lt_u, y=.data$Lls_u), color="red") + ggplot2::labs(x="Log-Age", y="Log of Log Survival")
        g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$Lt_l, y=.data$Lls_l, color="blue")) + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$Lt_l, y=.data$Lls_l), color="blue")
        g <- g + ggplot2::scale_color_manual(name = all_names[fir_KM], values = c("Above Mean" = "red", "Below Mean" = "blue")) 
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
#' @param fir term number for the intial term used for models of the form T0*f(Ti) in which the order matters
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
    for (fir_KM in 1:length(all_names)){
        lfir <- c(all_names[fir_KM])
        uniq <- unlist(unique(df[,lfir, with = FALSE]), use.names=FALSE)
        #
        e <- cox_ph_plot(Term_n, tform, a_n, er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type ,length(uniq))
        x <- e$x
        y <- e$y
        #
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
#            jpeg(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""), units="in", width=5, height=5, res=1200)
#            plot(x,y, type="l", xlab=all_names[fir_KM],ylab="Relative Risk",col='black',ylim=c(min(yl),max(yu)))
#            lines(xl,yl,col='black')
#            lines(xu,yu,col='black')
#            dev.off()
            #
            g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$x, y=.data$y)) + ggplot2::geom_line(color="black") + ggplot2::labs(x=all_names[fir_KM], y="Relative Risk")
            g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xl, y=.data$yl), color="black")
            g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xu, y=.data$yu), color="black")
            ggplot2::ggsave(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
            #
        } else {
#            jpeg(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""), units="in", width=5, height=5, res=1200)
#            plot(x,y, type="p", xlab=all_names[fir_KM],ylab="Relative Risk",col='black',ylim=c(min(yl),max(yu)))
#            lines(xl,yl,col='black',type="b")
#            lines(xu,yu,col='black',type="b")
#            dev.off()
            g <- ggplot2::ggplot(dft,ggplot2::aes(x=.data$x, y=.data$y)) + ggplot2::geom_point(color="black") + ggplot2::labs(x=all_names[fir_KM], y="Relative Risk")
            g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xl, y=.data$yl), color="black")
            g <- g + ggplot2::geom_line(data=dft, ggplot2::aes(x=.data$xu, y=.data$yu), color="black")
            ggplot2::ggsave(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""),device="jpeg",dpi="retina")
            #
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
#' @param fir term number for the intial term used for models of the form T0*f(Ti) in which the order matters
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
        if (verbose){
            print("plotting survival data")
#            print(t)
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
