#' Automatically assigns missing values in listed columns
#' \code{Replace_Missing} checks each column and fills in NA values
#'
#' @param df datatable being replaced in
#' @param name_list vector of string column names to check
#' @param MSV value to replace na with, same used for every column used
#' @param verbose verbosity for printing which columns needed fills
#'
#' @return returns a filled datatable
#' @export
#'
Replace_Missing <- function(df,name_list,MSV,verbose=FALSE){
    if (is.na(MSV)){
        if (verbose){
            print("The missing-value replacement is also NA")
        }
        stop()
    }
    for (j in name_list){
        #
        if (j %in% names(df)){
            ;
        } else {
            if (verbose){
                print(paste(j," missing from column names",sep=""))
            }
            stop()
        }
        #
        if (sum(is.na(df[[j]]))){
            set(df,which(is.na(df[[j]])),j,MSV)
            if (verbose){
                print(paste("Column ",j," had replaced values",sep=""))
            }
        }
    }
    return (df)
}


#' Automatically assigns missing control values
#' \code{Def_Control} checks and assigns default values
#'
#' @param control list of control parameters
#'
#' @return returns a filled list
#' @export
#'
Def_Control <- function(control){
    control_def=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1,"keep_strata"=FALSE,"Ncores"=detectCores())
    for (nm in names(control_def)){
        if (nm %in% names(control)){
            ;
        } else {
            control[nm] = control_def[nm]
        }
    }
    return (control)
}


#' Calculates Full Parameter list for Special Dose Formula
#' \code{Linked_Dose_Formula} Calculates all parameters for linear-quadratic and linear-exponential linked formulas
#'
#' @param tforms list of formula types
#' @param paras list of formula parameters
#' @param verbose verbosity argument, for error returns
#'
#' @return returns list of full parameters
#' @export
#'
Linked_Dose_Formula <- function(tforms,paras,verbose=FALSE){
    full_paras <- list()
    for (nm in names(tforms)){
        if (tforms[nm]=="quad"){
            plist <- unlist(paras[nm],use.names=FALSE)
            a0 <- plist[1]
            y <- plist[2]
            if (a0 < 0 ){
                if (verbose){
                    print("a0 arguement was negative")
                }
                stop()
            }
            if (is.numeric(a0)){
                ;
            } else {
                if (verbose){
                    print("a0 arguement was not a number")
                }
                stop()
            }
            if (is.numeric(y)){
                ;
            } else {
                if (verbose){
                    print("threshold arguement was not a number")
                }
                stop()
            }
            a1 <- a0/2/y
            b1 <- a0*y/2
            full_paras[[nm]] <- c(y,a0,a1,b1)
        } else if (tforms[nm]=="exp"){
            plist <- unlist(paras[nm],use.names=FALSE)
            a0 <- plist[1]
            y <- plist[2]
            b1 <-plist[3]
            if (a0 < 0 ){
                if (verbose){
                    print("a0 arguement was negative")
                }
                stop()
            }
            if (is.numeric(a0)){
                ;
            } else {
                if (verbose){
                    print("a0 arguement was not a number")
                }
                stop()
            }
            if (is.numeric(y)){
                ;
            } else {
                if (verbose){
                    print("threshold arguement was not a number")
                }
                stop()
            }
            if (is.numeric(b1)){
                ;
            } else {
                if (verbose){
                    print("exponential arguement was not a number")
                }
                stop()
            }
            c1 <- log(a0)-log(b1)+b1*y
            a1 <- a0*y+exp(c1-b1*y)
            full_paras[[nm]] <- c(y,a0,a1,b1,c1)
        }
    }
    return (full_paras)
}

#' Calculates The Additional Parameter For a linear-exponential formula with known maximum
#' \code{Linked_Lin_Exp_Para} Calculates what the additional parameter would be for a desired maximum
#'
#' @param y point formula switch
#' @param a0 linear slope
#' @param a1_goal exponential maximum desired
#' @param verbose verbosity argument, for error returns
#'
#' @return returns parameter used by Colossus
#' @export
#'
Linked_Lin_Exp_Para <- function(y,a0,a1_goal,verbose=FALSE){
    b1 <- 10
    lr <- 1.0
    if (a0 < 0 ){
        if (verbose){
            print("a0 arguement was negative")
        }
        stop()
    }
    if (a1_goal > y*a0){
        ;
    } else {
        if (verbose){
            print("goal is too low")
        }
        stop()
    }
    iter_i <- 0
    while (iter_i<100){
        iter_i = iter_i + 1
        c1=log(a0/b1)+b1*y
        a1=a0*y+exp(c1-b1*y)
        a_dif <- a1-a1_goal
#        if (abs(a_dif)<1e-3){
#            stop()
#        }
        if (abs(a_dif)<1e-3){
            break   
        }
        da1=-1/b1*exp(c1-b1*y)
        db1=(a1_goal-a1)/da1#-b1*a1_goal*np.exp(b1*y-c1)+b1*a0*y*np.exp(b1*y-c1)+b1
        if (-1*db1 > b1){
            db1=-0.9*b1
        }
        b1 = b1 + lr*db1
    }
    return (b1)
}



#' Opens a file for reading
#' \code{Open_File} uses user provided file name, time columns, and event column to load and sort.
#'
#' @param fname a string file name
#' @param time1 a column name defining the start of time periods, sorted third
#' @param time2 a column name defining the end of time periods, sorted first
#' @param event a column name defining the events, sorted second
#'
#' @return returns a sorted dataframe
#' @export
#'
Open_File <- function(fname,time1="age_start", time2="age_exit", event="cases"){
    colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=TRUE,fill=TRUE)
    setkeyv(df, c(time2, event,time1))
    return (df)
}

#' Splits a parameter into factors
#' \code{factorize} uses user provided list of columns to define new parameter for each unique value and update the data.table.
#' Not for interaction terms
#'
#' @param df a data.table containing the columns of interest
#' @param col_list an array of column names that should have factor terms defined
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#' @export
#'
factorize <-function(df,col_list){
    cols <- c()
    for (i in 1:length(col_list)){
        col <- col_list[i]
        x <- sort(unlist(as.list(unique(df[,col, with = FALSE])),use.names=FALSE))
    #    print(x)
        for (i in x){
            newcol <- c(paste(col,i,sep="_"))
            if (sum(df[,col, with = FALSE]==i)>0){
                df[, newcol] <- 1*(df[,col, with = FALSE]==i)
                cols <- c(cols, newcol)
            }
        }
    }
#    print(df)
    list('df'=df, 'cols'=cols)
}

#' Defines Interactions
#' \code{interact_them} uses user provided interactions define interaction terms and update the data.table. assumes interaction is "+" or "*" and applies basic anti-aliasing to avoid duplicates
#'
#' @param df a datatable containing the columns of interest
#' @param interactions array of strings, each one is of form term1?*?term2" for term1 interaction of type * with term2, "?" dlimits
#' @param new_names list of interaction names to use instead of default, default used if entry is ''
#' @param verbose verbosity
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#' @export
#'
interact_them <- function(df,interactions,new_names,verbose=FALSE){
    cols <- c()
    for (i in 1:length(interactions)){
        interac <- interactions[i]
        formula <- unlist(strsplit(interac,"[?]"),use.names=FALSE)
        if (length(formula)!=3){
            if (verbose){
                print(paste("Iteration:",interac,"has incorrect length of",length(formula),"but should be 3."))
            }
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
                df[, newcol] <- df[,col1, with = FALSE] + df[,col2, with = FALSE]
                cols <- c(cols, newcol)
            } else if (formula[2]=="*"){
                df[, newcol] <- df[,col1, with = FALSE] * df[,col2, with = FALSE]
                cols <- c(cols, newcol)
            } else {
                if (verbose){
                    print(paste("Incorrect operation of",formula[2]))
                }
                stop()
            }
        }
    }
    list('df'=df, 'cols'=cols)
}

#' Defines the likelihood ratio test
#' \code{Likelihood_Ratio_Test} uses two models and calculates the ratio
#'
#' @param alternative_model the new model of interest
#' @param null_model a model to compare against
#'
#' @return returns the score statistic
#' @export
#'
Likelihood_Ratio_Test <- function(alternative_model, null_model){
    if (("LogLik" %in% names(alternative_model))&&("LogLik" %in% names(null_model))){
        return (2*(unlist(alternative_model["LogLik"],use.names=FALSE) - unlist(null_model["LogLik"],use.names=FALSE)))
    }
    stop()
    return (NULL)
}


#' Calculates a chi-squared value, used in a non-functioning bounds formula
#' \code{get_conf_int} uses an interval parameter to return the score statistic
#'
#' @param alpha decimal 1 - (confidence interval)
#'
#' @return returns the score statistic
#' @export
#'
get_conf_int <-function(alpha=0.95){
    q1 <- qchisq(1-alpha, df=1)
    return (q1)
}


#' checks for duplicated column names
#' \code{Check_Dupe_Columns} checks for duplicated columns, columns with the same values, and columns with 1 value. Currently not updated for multi-terms
#'
#' @param df dataframe of data to use as reference
#' @param cols columns to check
#' @param term_n term numbers
#' @param verbose verbosity
#'
#' @return returns the usable columns
#' @export
#'
Check_Dupe_Columns <- function(df,cols,term_n,verbose=FALSE){
    ##
    if (length(cols)>1){
        features_pair <- combn(cols, 2, simplify = F) # list all column pairs
        terms_pair <- combn(term_n, 2, simplify = F) # list all term pairs
        toRemove <- c() # init a vector to store duplicates
        for(pair_n in 1:length(features_pair)) { # put the pairs for testing into temp objects
            pair <- unlist(features_pair[pair_n])
            term <- unlist(terms_pair[pair_n])
            f1 <- pair[1]
            f2 <- pair[2]
            t1 <- term[1]
            t2 <- term[2]
            if (t1==t2){
                df[,get(f1)]
                df[,get(f2)]
                if (!(f1 %in% toRemove) & !(f2 %in% toRemove)) {
                    if (all(df[[f1]] == df[[f2]])) { # test for duplicates
                        if (verbose){
                            print(paste(f1, " and ", f2, " are equals.",sep=""))
                        }
                        toRemove <- c(toRemove, f2) # build the list of duplicates
                    }
                    if (min(df[[f2]])==max(df[[f2]])){
                        if (min(df[[f2]])==0){
                            toRemove <- c(toRemove, f2) # remove zero values
                        }
                    }
                }
            }
        }
        newcol <- setdiff(cols, toRemove)
        if (length(newcol)==1){
            if (min(df[,newcol, with = FALSE])==max(df[,newcol, with = FALSE])){
                return(c())
            } else {
                return(newcol)
            }
        }
        return(newcol)
    } else if (length(cols)==1){
        if (min(df[,cols, with = FALSE])==max(df[,cols, with = FALSE])){
            return(c())
        } else {
            return(cols)
        }
    } else {
        return(c())
    }
    return(c())
}

#' Applies time duration truncation limits
#' \code{Check_Trunc} creates columns to use for truncation
#'
#' @param df dataframe of data to use as reference
#' @param ce columns to check for truncation
#' @param verbose verbosity argument, for error returns
#'
#' @return returns the updated data and time period columns
#' @export
#'
Check_Trunc <- function(df,ce,verbose=FALSE){
    if (ce[1]=="%trunc%"){
        if (ce[2]=="%trunc%"){
            if (verbose){
                print("Both endpoints are truncated, not acceptable")
            }
            stop()
        }
        tname <- ce[2]
        tmin <- min(df[,get(tname)])-1
        if (!("right_trunc" %in% names(df))){
            df[,':='(right_trunc=tmin)]
        }
        ce[1]="right_trunc"
    } else if (ce[2]=="%trunc%") {
        tname <- ce[1]
        tmax <- max(df[,get(tname)])+1
        if (!("left_trunc" %in% names(df))){
            df[,':='(left_trunc=tmax)]
        }
        ce[2]="left_trunc"
    }
    return (list('df'=df,'ce'=ce))
}

#' Applies time depedence to parameters
#' \code{gen_time_dep} generates a new dataframe with time dependent covariates by applying a grid in time
#'
#' @param df dataframe of data to use as reference
#' @param time0 starting time column
#' @param time1 ending time column
#' @param event event column
#' @param iscox boolean if rows not at event times should be kept
#' @param dt spacing in time for new rows
#' @param new_names vector of new column names for the time dependent columns
#' @param dep_cols columns that are not needed in the new dataframe
#' @param func_form vector of functions to apply to each time-dependent covariate. Of the form func(df, time) returning a vector of the new column value
#' @param fname filename used for new dataframe
#' @param tform list of string function identifiers, used for linear/step
#'
#' @return returns the updated dataframe
#' @export
#'
gen_time_dep <- function(df, time0, time1, event, iscox, dt, new_names, dep_cols, func_form,fname, tform){
    dfn <- names(df)
    ce <- c(time0,time1,event)
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    time0 <- ce[1]
    time1 <- ce[2]
    dfn_same <- dfn[!(dfn %in% dep_cols)]
    dfn_dep <- c()
    for (i in 1:length(new_names)){
        name0 <- paste(new_names[i],0,sep="_")
        name1 <- paste(new_names[i],1,sep="_")
        func <- func_form[i]
        df[, name0] = lapply(func, function(f) f(df, time0))
        df[, name1] = lapply(func, function(f) f(df, time1))
        dfn_dep <- c(dfn_dep, name0, name1)
    }
    #
    #
    dfn_time <- c(time0, time1)
    dfn_event <- c(event)
    dfn_same <- dfn_same[!(dfn_same %in% dfn_time)]
    dfn_same <- dfn_same[!(dfn_same %in% dfn_event)]
    ##
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time1, with = FALSE]), use.names=FALSE))
    if (iscox){
        #
        df <- df[get(time1)>=min(tu),]
        df <- df[get(time0)<=max(tu),]
        #
    }
    ##
    x_time = as.matrix(df[,dfn_time, with = FALSE])
    x_dep = as.matrix(df[,dfn_dep, with = FALSE])
    x_same = as.matrix(df[,dfn_same, with = FALSE])
    x_event = as.matrix(df[,dfn_event, with = FALSE])
    #
    if (grepl(".csv", fname, fixed = TRUE)){
        ;
    } else {
        fname <- paste(fname,".csv",sep="_")
    }
    Write_Time_Dep(x_time, x_dep, x_same, x_event, dt, fname,tform,tu,iscox)
    df_new <- fread(fname,data.table=TRUE,header=FALSE,col.names=c(time0,time1,new_names,dfn_same,event))
    setkeyv(df_new, c(time1, event))
    return (df_new)
}

#' Automates creating a date difference column
#' \code{Date_Shift} generates a new dataframe with a column containing time difference in a given unit
#'
#' @param df dataframe of data to use as reference
#' @param dcol0 list of starting month, day, and year
#' @param dcol1 list of ending month, day, and year
#' @param col_name new column name
#' @param units time unit to use
#'
#' @return returns the updated dataframe
#' @export
#'
Date_Shift <- function(df, dcol0, dcol1, col_name, units="days"){
    def_cols <- names(df)
    #
    df$dt0 <- paste(df[[match(dcol0[1],names(df))]],df[[match(dcol0[2],names(df))]],df[[match(dcol0[3],names(df))]],sep="-")
    df$dt1 <- paste(df[[match(dcol1[1],names(df))]],df[[match(dcol1[2],names(df))]],df[[match(dcol1[3],names(df))]],sep="-")
    #
    # TO NOT ENCOUNTER DAYLIGHT SAVINGS ISSUES, THE UTC TIMEZONE IS USED
    # IF NOT USED THEN RESULTS MAY HAVE TIMES OFF BY 1/24 decimals
    df[, col_name] = difftime(strptime(df$dt1, format = "%m-%d-%Y",tz = 'UTC'), strptime(df$dt0,  format = "%m-%d-%Y"), units = units,tz = 'UTC')
    def_cols <- c(def_cols, col_name)
    return (df[,def_cols,with=FALSE])
}

#' Automates creating a date since a reference column
#' \code{Time_Since} generates a new dataframe with a column containing time since a reference in a given unit
#'
#' @param df dataframe of data to use as reference
#' @param dcol0 list of ending month, day, and year
#' @param tref reference time in date format
#' @param col_name new column name
#' @param units time unit to use
#'
#' @return returns the updated dataframe
#' @export
#'
Time_Since <- function(df, dcol0, tref, col_name, units="days"){
    def_cols <- names(df)
    #
    df$dt0 <- paste(df[[match(dcol0[1],names(df))]],df[[match(dcol0[2],names(df))]],df[[match(dcol0[3],names(df))]],sep="-")
    #
    #
    df[, col_name] = lapply(df$dt0, function(x) (difftime(strptime(x,  format = "%m-%d-%Y",tz = 'UTC'), tref, units = units)))
    def_cols <- c(def_cols, col_name)
    return (df[,def_cols,with=FALSE])
}
