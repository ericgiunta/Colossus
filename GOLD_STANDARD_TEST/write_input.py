script_temp='''library(Colossus)
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

path <- "{base}_LL_{extension}.parquet"
path <- path.expand(path)
path <- normalizePath(path)
xdf <- pandas$read_parquet(path)
df0 <- as.data.table(xdf)


setkeyv(df,  c("seq_id"))
setkeyv(df0, c("seq_id"))

df$lung     <- df0$event
#df$age_exit     <- df0$exit_adjust

df <- df %>% mutate_at(vars(age_entry, age_exit), funs(compx_round(., 2)))

df <- df[lung<2]

time1 <- "age_entry"
time2 <- "age_exit"
event <- "lung"


col_list <- c("sex",'Urban','Work_Group','Smoking')

temp <- factorize(df,col_list)
df <- temp$df


names <- c('sex_1', 'Urban_1', 'Work_Group_1', 'Work_Group_2', 'Work_Group_3', 'Work_Group_4', 'Smoking_1', 'dose')

Term_n <- rep(0,length(names))
tform <- rep("loglin",length(names))
keep_constant <- rep(0,length(names))
a_n <- rep(0.1,length(names))
modelform <- "M"
fir <- 0
der_iden <- 0
control=list('lr' = 0.75,'maxiter' = 50,'halfmax' = 5,'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
print(paste("df404 ",nrow(df),sep=""))
print(paste("df405 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
e <- RunCoxRegression_Basic(df, time1, time2, event, names, keep_constant, a_n, der_iden, control)
#
print(paste("df406 ",difftime(Sys.time(),stime,unit="secs"),sep=""))
print("output: {base}_{extension}")
print(e)
'''

bash_temp0='''#!/bin/bash -l

## Specify RAM needed per core.  Default is 1G.
#SBATCH --mem=70G   # Memory per core, use --mem= for memory per node
## Specify maximum runtime.  Default is 1 hour (1:00:00)
#SBATCH --time=00-0:30:00   # Use the form DD-HH:MM:SS
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

bash_temp1='''#!/bin/bash -l

## Specify RAM needed per core.  Default is 1G.
#SBATCH --mem=200G   # Memory per core, use --mem= for memory per node
## Specify maximum runtime.  Default is 1 hour (1:00:00)
#SBATCH --time=00-1:30:00   # Use the form DD-HH:MM:SS
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

def write_gen(base,extension):
    fname="{}_{}.r".format(base,extension)
    f=open(fname,'w')
    f.write(script_temp.replace("{base}",base).replace("{extension}",extension))
    f.close()
    
    fname="{}_{}.sh".format(base,extension)
    f=open(fname,'w')
    allowed=[]
    for i in ["5050","2575","7525"]:
        for j in ["44051","44052","44053","44054"]:
            allowed.append(i+"_"+j)
    if ("{base}_{extension}".replace("{base}",base).replace("{extension}",extension) in allowed):
        f.write(bash_temp1.replace("{base}",base).replace("{extension}",extension))
    else:
        f.write(bash_temp0.replace("{base}",base).replace("{extension}",extension))
    f.close()

def add_to_submit(base,extension,submit_extra):
    allowed=[]
    for i in ["5050","2575","7525"]:
        for j in ["44051","44052","44053","44054"]:
            allowed.append(i+"_"+j)
    if ("{base}_{extension}".replace("{base}",base).replace("{extension}",extension) in allowed):
        submit_extra+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
    return submit_extra

submit_test=''
submit_rest=''
submit_extra=''
for base in ["5050","2575","7525"]:
    for i in [1,2,3,4,5]:
        extension="{}4055".format(i)
        write_gen(base,extension)
        submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
        submit_extra = add_to_submit(base,extension,submit_extra)
    for i in [1,2,3,5]:
        extension="4{}055".format(i)
        write_gen(base,extension)
        submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
        submit_extra = add_to_submit(base,extension,submit_extra)
    for i in [1,2]:
        extension="44{}55".format(i)
        write_gen(base,extension)
        submit_test+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
        submit_extra = add_to_submit(base,extension,submit_extra)
    for i in [1,2,3,4]:
        extension="440{}5".format(i)
        write_gen(base,extension)
        submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
        submit_extra = add_to_submit(base,extension,submit_extra)
    for i in [1,2,3,4]:
        extension="4405{}".format(i)
        write_gen(base,extension)
        submit_rest+="\nsbatch {base}_{extension}.sh\nsleep 1\n".replace("{base}",base).replace("{extension}",extension)
        submit_extra = add_to_submit(base,extension,submit_extra)
    
f=open("submit_test.sh",'w')
f.write(submit_test)
f.close()

f=open("submit_rest.sh",'w')
f.write(submit_rest)
f.close()

f=open("submit_extra.sh",'w')
f.write(submit_extra)
f.close()









