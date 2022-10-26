#' Performs basic Cox Proportional Hazards regression
#' \code{RunCoxRegression} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
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
#' @param fir term number for the intial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#'
RunCoxRegression <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    dfend <- df[get(event)==1, ]
    tu <- unlist(unique(dfend[,time2, with = FALSE]))
    if (length(tu)==0){
        stop()
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
    a_n0 <- copy(a_n)
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    a_n <- a_n0
    return (e)
}

#' Calculates proportional hazard for a reference vector
#' \code{RunCoxRegression} uses user provided data, time/event columns, vectors specifying the model, and options to calculate risk for a reference
#'
#' @param df data used to calculate PH
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the intial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters for verbosity and tie method
#'
#' @return returns a list of the final results
#' @export
#'
Cox_Relative_Risk <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control){
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    #
    #
    e <- cox_ph_risk_sub(Term_n, tform, a_n, dfc, x_all,  fir, modelform, control, term_tot)
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with the null model
#' \code{RunCoxRegression} uses user provided data and time/event columns to calculate the log-likelihood with constant hazard ratio
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#'
RunCoxNull <- function(df, time1, time2, event,control){
    dfend <- df[get(event)==1, ]
    tu <- unlist(unique(dfend[,time2, with = FALSE]))
    if (length(tu)==0){
        stop()
    }
    ce <- c(time1,time2,event)
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    e <- cox_ph_null( control, as.matrix(df[,ce, with = FALSE]), tu)
    return (e)

}


#' Performs Cox Proportional Hazard model plots
#' \code{RunCoxPlots} uses user provided data, time/event columns, vectors specifying the model, and options to choose plots and saves plots
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
#' @param fir term number for the intial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param plot_options list of parameters controlling the plot options
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
#'
RunCoxPlots <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, plot_options){
    IDS <- base <- res  <- NULL
    verbose <- copy(plot_options$verbose)
    verbosec <- copy(control$verbose)
    dfend <- df[get(event)==1, ]
    tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
    if (length(tu)==0){
        stop()
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
    time1 <- ce[1]
    time2 <- ce[2]
    Plot_Type <- plot_options$type
    #
    control$verbose <- verbose
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    control$verbose <- verbosec
    b <- e$beta_0
    er <- e$Standard_Deviation
    #
    #
    # IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, vector<string> Plot_Type ,int uniq_v
    e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
    #
    t <- c()
    h <- c()
    ch <- c()
    surv <- c()
    dt <- 1
    if (verbose){
        print("writing survival data")
    }
    dft=data.table("time"=tu,"base"=e$baseline)
    for (i in tu){
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
    #
    if (Plot_Type[1]=="SURV"){
        if (plot_options$Martingale==TRUE){
            #
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
            dname <- plot_options$dose_col
            e_i <- df[,get(event)]
            dfr=data.table("Risks"=e$Risks,"ch_e"=ch_e,"ch_s"=ch_s,"e"=e_i,"IDS"=df$studyID,"time"=time_e,"doses"=df[, dname, with = FALSE])
            dfr$res = dfr$e - (dfr$ch_e-dfr$ch_s) * dfr$Risks
            #
            Martingale_Error <- dfr[, lapply(.SD,sum), by=IDS]
            times <- dfr[, lapply(.SD,max), by=IDS]
            ##
            jpeg(paste('dose_sum_by_max',Plot_Type[2],'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
            plot(times[,doses],Martingale_Error[,doses], xlab="Max Dose",ylab="Sum of Doses",type='p')
            dev.off()
            ##
            jpeg(paste('martin_plot',Plot_Type[2],'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
            smoothScatter(times[,time],Martingale_Error[,res], xlab="age",ylab="Martingale Residuals",nbin=100, colramp= colorRampPalette(c("white","red")))
            dev.off()
            jpeg(paste('martin_plot_dose',Plot_Type[2],'.jpg',sep="_"), units="in", width=5, height=5, res=1200)
            smoothScatter(times[doses>0,doses],Martingale_Error[doses>0,res], xlab="Dose",ylab="Martingale Residuals",nbin=200, colramp= colorRampPalette(c("white","red")))
            dev.off()
            #
        }
        if (plot_options$surv_curv==TRUE){
            #
            if (verbose){
                print("Survival Plots")
            }
            #
            jpeg(paste('ch_plot',Plot_Type[2],'.jpg',sep="_"), units="in", width=5, height=5, res=1200)
            plot(t,ch, type="s", xlab="age",ylab="Cumulative Hazard")
            dev.off()
            #
            jpeg(paste('surv_plot',Plot_Type[2],'.jpg',sep="_"), units="in", width=5, height=5, res=1200)
            plot(t,surv, type="s", xlab="age",ylab="Survival")
            dev.off()
            #
            for (i in 1:length(a_n)){
                a_n[i] <- b[i] - er[i]
            }
            e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
            dfl=data.table("time"=tu,"base"=e$baseline)
            #
            for (i in 1:length(a_n)){
                a_n[i] <- b[i] + er[i]
            }
            e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
            dfu=data.table("time"=tu,"base"=e$baseline)
            #
#            hl <- c()
#            chl <- c()
#            hu <- c()
#            chu <- c()
#            if (verbose){
#                print("writing survival data")
#            }
#            dft=data.table("time"=tu,"base"=e$baseline)
#            for (i in tu){
#                temp <- sum(dfl[time<i, base])
#                if (length(hl)==0){
#                    hl <- c(temp)
#                } else {
#                    hl <- c(hl, temp-chl[length(chl)])
#                }
#                chl <- c(chl, temp)
#                temp <- sum(dfu[time<i, base])
#                if (length(hu)==0){
#                    hu <- c(temp)
#                } else {
#                    hu <- c(hu, temp-chu[length(chu)])
#                }
#                chu <- c(chu, temp)
#            }
            jpeg(paste("H_plot",Plot_Type[2],".jpg",sep="_"), units="in", width=5, height=5, res=1200)
            plot(t,h, xlab="age",ylab="Hazard",col='black',type='s')
#            lines(t,hl,col='blue',type='p')
#            lines(t,hu,col='red',type='p')
#            legend("topright", legend=c("Estimate", "Low Estimate","High Estimate"), col=c("blue","black", "red"), lty=1:2, cex=0.8)
            dev.off()
        }
        if (plot_options$KM==TRUE){
            #
            if (verbose){
                print("KM Plots")
            }
            studyID <- plot_options$studyID
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
                jpeg(paste("KM_",fir_KM,"_",Plot_Type[2],".jpg",sep=""), units="in", width=5, height=5, res=1200)
                plot(t_u,n_u, xlab=all_names[fir_KM],ylab="Survival Fraction",ylim=c(min(n_l),max(n_u)),col='red',type='l')
                lines(t_l,n_l,col='blue')
                legend("topright", legend=c("Below Mean", "Above Mean"), col=c("red", "blue"), lty=1:2, cex=0.8,title=all_names[fir_KM])
                dev.off()
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
                control$verbose <- verbose
                e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_l[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
                control$verbose <- verbosec
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
                LLs_l <- log(-Ls)
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
                control$verbose = verbose
                e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_u[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
                control$verbose = verbosec
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
                LLs_u <- log(-Ls)
                Lt_u <- log(t)
                ##
                #
                if (verbose){
                    print("Survival Plots")
                }
                #
                jpeg(paste('log_log_surv_plot',fir_KM,"_",Plot_Type[2],'.jpg',sep='_'), units="in", width=5, height=5, res=1200)
                plot(Lt_l,LLs_l, xlab="Log-Age",ylab="Log of Log Survival",col='red',type='l')
                lines(Lt_u,LLs_u,col='blue')
                legend("topleft", legend=c("Below Mean", "Above Mean"), col=c("red", "blue"), lty=1:2, cex=0.8,title=all_names[fir_KM])
                dev.off()
            }
            #
        }
    } else if (Plot_Type[1]=="RISK"){
        fir_KM=0
        dfend <- df[get(event)==1, ]
        #
        tu <- unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE)
        print("start risk calculations")
        for (fir_KM in 1:length(all_names)){
            lfir <- c(all_names[fir_KM])
            uniq <- unlist(unique(df[,lfir, with = FALSE]), use.names=FALSE)
            #
            e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df_l[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type ,length(uniq))
            x <- e$x
            y <- e$y
            #
            a_n[fir_KM] <- b[fir_KM] - er[fir_KM]
            #
            e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , length(uniq))
            xl <- e$x
            yl <- e$y
            #
            a_n[fir_KM] <- b[fir_KM] + er[fir_KM]
            #
            e <- cox_ph_plot(Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , length(uniq))
            xu <- e$x
            yu <- e$y
            if (length(uniq)>100){
                jpeg(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""), units="in", width=5, height=5, res=1200)
                plot(x,y, type="l", xlab=all_names[fir_KM],ylab="Relative Risk",col='black',ylim=c(min(yl),max(yu)))
                lines(xl,yl,col='black')
                lines(xu,yu,col='black')
                dev.off()
                #
            } else {
                jpeg(paste("risk_plot_",fir_KM,"_",Plot_Type[2],".jpg",sep=""), units="in", width=5, height=5, res=1200)
                plot(x,y, type="p", xlab=all_names[fir_KM],ylab="Relative Risk",col='black',ylim=c(min(yl),max(yu)))
                lines(xl,yl,col='black',type="b")
                lines(xu,yu,col='black',type="b")
                dev.off()
                #
            }
            #
        }
    }
    return ("Passed")
}






















