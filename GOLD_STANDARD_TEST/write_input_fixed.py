script_temp_EXP='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','CONST','CONST','dose','dose')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'step_slope','step_int','loglin_slope',"loglin_top")
keep_constant <- rep(0,length(names))
keep_constant[1:10] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6, {const_val0}, 0, {const_val1}, 0.001)#rep(0.1,length(names))
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,9),1,0),'rmin'=c(0,0,0,0,0,0,0,0,0,0,-4.5/30),'rmax'=c(0,0,0,0,0,0,0,0,0,0,1.6/30))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''


script_temp_L_EXP='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df

df$dose2 <- df$dose * df$dose

names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','dose','dose2')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'plin','loglin')
keep_constant <- rep(0,length(names))
keep_constant[1:7] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6,0.0, -0.001)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,7),0,0),'rmin'=c(0,0,0,0,0,0,0,-1/30,-0.005),'rmax'=c(0,0,0,0,0,0,0,1/30,0.002))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

script_temp_L_Q='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','CONST','CONST','dose','dose','dose')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'step_slope','step_int',"lin_slope",'lin_int','quad_slope')
keep_constant <- rep(0,length(names))
keep_constant[8] <- 1
keep_constant[1:9] <- 1
keep_constant[11] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6, 1, 0 , 0.0,0.0,0.0)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,9),0,0,0),'rmin'=c(0,0,0,0,0,0,0,0,0,-1/30,0.0,-1/900),'rmax'=c(0,0,0,0,0,0,0,0,0,1/30,0.0,1/900))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

script_temp_Q='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','CONST','CONST','dose')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'step_slope','step_int','quad_slope')
keep_constant <- rep(0,length(names))
keep_constant[8] <- 1
keep_constant[1:9] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6, {const_val0}, 0, 0.001)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,9),1,0),'rmin'=c(0,0,0,0,0,0,0,0,0,-1/900),'rmax'=c(0,0,0,0,0,0,0,0,0,5/900))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

script_temp_ST='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','CONST','CONST','dose','dose')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'step_slope','step_int','step_slope','step_int')
keep_constant <- rep(0,length(names))
keep_constant[8] <- 1
keep_constant[1:9] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6, {const_val0}, 0, 0.1,10)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,9),0,0),'rmin'=c(0,0,0,0,0,0,0,0,0,-1.0,-10),'rmax'=c(0,0,0,0,0,0,0,0,0,1.0,20))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

script_temp_LNT='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','CONST','CONST','dose','dose')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'step_slope','step_int',"lin_slope","lin_int")
keep_constant <- rep(0,length(names))
keep_constant[8] <- 1
keep_constant[1:9] <- 1
keep_constant[11] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6, {const_val0}, 0, 0.0,0)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,9),0,0),'rmin'=c(0,0,0,0,0,0,0,0,0,-1/30,0.0),'rmax'=c(0,0,0,0,0,0,0,0,0,1/30,0.0))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

script_temp_LT='''library(Colossus)
library(data.table)
library(reticulate)
pandas <- import("pandas")
library(dplyr)
library(ggplot2)

stime <- Sys.time()

compx_round <- function(x,digits){
    a <- floor(x)
    b <- x-a
    b <- b*10^digits
    b <- ceiling(b)
    b <- b*10^(-1*digits)
    a+b
}

path <- "{base}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df <- as.data.table(xdf)

path <- "{base}_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 1)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1','CONST','CONST','dose','dose')

Term_n <- rep(0,length(names))
tform <- c(rep("loglin",7),'step_slope','step_int',"lin_slope","lin_int")
keep_constant <- rep(0,length(names))
keep_constant[8] <- 1
keep_constant[1:9] <- 1
a_n <- c(0.1,0.1,0.69,0.1,-0.1,-0.69,0.6, {const_val0}, 0, 0.0,10)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.95,'maxiter' = 10,'halfmax' = 3,'epsilon' = 1e-5,'dbeta_max' = 0.95,'deriv_epsilon' = 1e-5, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=2.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("maxiter"=2, 'verbose'=TRUE ,"guesses"=10,'strata'=FALSE, 'guess_constant'=c(rep(1,9),0,0),'rmin'=c(0,0,0,0,0,0,0,0,0,-1/10,-10),'rmax'=c(0,0,0,0,0,0,0,0,0,1/10,20))
Strat_Col='sex_1'
print(paste("df401 ",nrow(df),sep=""))
e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

bash_temp1='''#!/bin/bash -l

## Specify RAM needed per core.  Default is 1G.
#SBATCH --mem=450G   # Memory per core, use --mem= for memory per node
## Specify maximum runtime.  Default is 1 hour (1:00:00)
#SBATCH --time=00-3:00:00   # Use the form DD-HH:MM:SS
## Send email when job is aborted (BEGIN), begins (FAIL), or ends (END) (ALL does all three)
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END
## Email address
#SBATCH --mail-user=egiunta@ksu.edu
## Name the job
#SBATCH --job-name={base}_{extension}
## Number of nodes you'd like to target for this run
#SBATCH --nodes=1
## Number of parallel operations taking place at the same time on the same node (number of cores utilized per node)
#SBATCH --ntasks-per-node=32
#SBATCH --gres=killable:1
#SBATCH --partition=ksu-mne-bahadori.q

module load R/4.2.1-foss-2022a
source ~/virtualenvs/ret_virtual/bin/activate
export RETICULATE_PYTHON=/homes/egiunta/virtualenvs/ret_virtual/bin/python

Rscript {base}_{extension}.r &> {base}_{extension}.txt
'''

def write_gen(base,extension,const_val0,const_val1):
    const_val0 = str(const_val0)
    const_val1 = str(const_val1)
    fname="{}_{}.r".format(base,extension)
    f=open(fname,'w')
    if "L_EXP" in extension:
        f.write(script_temp_L_EXP.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0))
    elif "EXP" in extension:
        f.write(script_temp_EXP.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0).replace("{const_val1}",const_val1))
    elif "L_Q" in extension:
        f.write(script_temp_L_Q.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0))
    elif "Q" in extension:
        f.write(script_temp_Q.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0))
    elif "LNT" in extension:
        f.write(script_temp_LNT.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0))
    elif "LT" in extension:
        f.write(script_temp_LT.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0))
    elif "ST" in extension:
        f.write(script_temp_ST.replace("{base}",base).replace("{extension}",extension).replace("{const_val0}",const_val0))
    f.close()
    
    fname="{}_{}.sh".format(base,extension)
    f=open(fname,'w')
    f.write(bash_temp1.replace("{base}",base).replace("{extension}",extension))
    f.close()
    return


submit_test=''
submit_rest=''
#submit_extra=''
for base in ["5050","2575","7525"]:
    # for i in [1,2,4,5]:
    #     const_val0=1
    #     const_val1=0
    #     extension="L_EXP_4405{}".format(i)
    #     write_gen(base,extension,const_val0,const_val1)
    #     if (i==1):
    #         submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    #     else:
    #         submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
#    for i in [1,2,4,5]:
#        const_val0=1
#        const_val1=0
#        extension="LNT_4405{}".format(i)
#        write_gen(base,extension,const_val0,const_val1)
#        if (i==1):
#            submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
#        else:
#            submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    for i in [1,2,4,5]:
       const_val0=1
       const_val1=0
       extension="LT_4405{}".format(i)
       write_gen(base,extension,const_val0,const_val1)
       if (i==1):
           submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
       else:
           submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    # for i in [1,2,4,5]:
    #     const_val0=1
    #     const_val1=0
    #     extension="L_Q_4405{}".format(i)
    #     write_gen(base,extension,const_val0,const_val1)
    #     if (i==1):
    #         submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    #     else:
    #         submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    for i in [1,2,4,5]:
       const_val0=1
       const_val1=0
       extension="ST_4405{}".format(i)
       write_gen(base,extension,const_val0,const_val1)
       if (i==1):
           submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
       else:
           submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    # for i in [1,2,4,5]:
    #     if i in [1,2]:
    #         const_val0=2
    #         const_val1=-1
    #     elif i in [4,5]:
    #         const_val0=0
    #         const_val1=1
    #     extension="EXP_4405{}".format(i)
    #     dose_tforms=""
    #     write_gen(base,extension,const_val0,const_val1)
    #     if (i==1):
    #         submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    #     else:
    #         submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
#    for i in [1,2,4,5]:
#        const_val0=1
#        const_val1=0
#        extension="Q_4405{}".format(i)
#        write_gen(base,extension,const_val0,const_val1)
#        if (i==1):
#            submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
#        else:
#            submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)

    
f=open("submit_test.sh",'w')
f.write(submit_test)
f.close()

f=open("submit_rest.sh",'w')
f.write(submit_rest)
f.close()

#f=open("submit_extra.sh",'w')
#f.write(submit_extra)
#f.close()









