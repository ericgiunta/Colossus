library(Colossus)
library(data.table)
library(survival)


data(cancer, package="survival")

df <- cancer
df$UserID <- 1:nrow(df)

df$status <- df$status - 1
df$sex <- df$sex - 1

t0 <- "%trunc%"
t1 <- "time"
event <- "status"

names <- c('age',"sex")
tform <- c("loglin", "loglin")
control <- list("Ncores"=2, 'maxiter'=50, 'verbose'=3)

#e <- RunCoxRegression_Basic(df, t0, t1, event, names, control = control)

a_n <- c(0.01701289, -0.51256478)
term_n <- c(0,0)
keep_constant <- c(0,0)
modelform <- 'M'
fir <- 0

plot_options <- list("type"=c("surv",paste(tempfile(),"run",sep="")), "studyid"="UserID", 'verbose'=4, "surv_curv"=T, 'martingale'=F,'strat_haz'=F,'km'=F)

#e <- RunCoxPlots(df, t0, t1, event, names, term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)

plot_options <- list("type"=c("surv",paste(tempfile(),"run",sep="")), "studyid"="UserID", 'verbose'=4, "surv_curv"=F, 'martingale'=T,'strat_haz'=F,'km'=F,'cov_cols'=c('age','sex'))
e <- RunCoxPlots(df, t0, t1, event, names, term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)

print(e)
plot_options <- list("type"=c("surv",paste(tempfile(),"run",sep="")), "studyid"="UserID", 'verbose'=4, "surv_curv"=F, 'martingale'=F,'strat_haz'=F,'km'=T)

plot_options <- list("type"=c("risk",paste(tempfile(),"run",sep="")), "studyid"="UserID", 'verbose'=4)
plot_options <- list("type"=c("schoenfeld",paste(tempfile(),"run",sep="")), "studyid"="UserID", 'verbose'=4)
