library(data.table)
library(parallel)
library(Colossus)

fname <- 'base_example.csv'
df <- fread(fname)
df$pyr <- df$exit - df$entry
# 
time1 <- "entry"
time2 <- "exit"
pyr <- "pyr"
event <- "event"
names <- c("dose0","dose1", "exit")
term_n <- c(0,0,0)
tform <- c("loglin","loglin","loglin")
keep_constant <- c(0,0,0)
#a_n <- c(0.2462, 5.020, -0.5909)
a_n <- c(-0.6067272,  5.0193899, -4.6564575)
modelform <- "M"
fir <- 0
der_iden <- 0
#
model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=FALSE, 'alpha'=0.1)
control <- list("ncores"=2,'lr' = 0.75,'maxiters' = c(50,50),'halfmax' = 5,'epsilon' = 0,'dbeta_max' = 0.5,'deriv_epsilon' =0, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=4, 'ties'='breslow','double_step'=1)
# 
e <- RunPoissonRegression_Omnibus(df,pyr, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
print(e)
stop()
sink(file="out.txt")
alpha_list <- c(0.75, 0.5, 1-0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
#alpha_list <- c( 0.005)
#v_lower <- c(0.10224367838276, -0.0366388051905475, -0.219114008730643, -0.344171753035397, -1.95449016956737, -2.02344865739314, -2.08312767847433, -2.1518391834648, -2.19804353276631)
#v_upper <- c(0.369809870835033, 0.458228372102368, 1.07214890234804, 0.575303932702712, 0.688777681184681, 0.756984966300005, 0.81554103308482, 0.882549146565759, 0.927406355126396)
#for (alpha_i in 1:length(alpha_list)){
#    alpha <- alpha_list[alpha_i]
#    a_n <- c(0.2462, 5.020,-0.599)
#    model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=0, 'manual'=FALSE)
#    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
#    a <- e$Parameter_Limits
#    print(paste(a[1], a[2], sep=" "))
#     
#}
# 
# v_lower <- c(-0.643365949558998, -0.677336655540846, -0.706075250211414, -0.718165409196492, -0.753647332793819, -0.773208334303991, -0.789018704115451, -0.806061085000755, -0.816875114954096)
# v_upper <- c(-0.521472203247917, -0.444964438813732, -0.327862977142017, -0.235044092073815, 2.91573713669059, 3.21014641617297, 3.48490803194128, 3.82648584413642, 4.07272009904963)
# for (alpha_i in 1:length(alpha_list)){
#   alpha <- alpha_list[alpha_i]
#   a_n <- c(0.2462, 5.020,-0.599)
#   model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=2, 'manual'=FALSE)
#   e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
#   a <- e$Parameter_Limits
#   print(paste(a[1], v_lower[alpha_i], a[2], v_upper[alpha_i], sep=" "))
#   
# }
# sink(NULL)
# stop()
# v_lower <- c(4.97252283680023, 4.93499451084736, 4.8980471564778, 4.88084912217379, 4.89656496320952, 4.8436813543801, 4.802128281116, 4.75611040136942, 4.72560413563056)
# v_upper <- c(5.06762896556864, 5.10561529673742, 5.14327069553938, 5.16088604907152, 5.21982880779017, 5.25778604710117, 5.29187760189931, 5.33258757876175, 5.36084782860162)
# for (alpha_i in 1:length(alpha_list)){
#     alpha <- alpha_list[alpha_i]
#     a_n <- c(0.2462, 5.020,-0.599)
#     model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=1, 'manual'=FALSE)
#     e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
#     a <- e$Parameter_Limits
#     print(paste(a[1] , v_lower[alpha_i], a[2] , v_upper[alpha_i], sep=" "))
# }
# # #
# v_lower <- c(-0.643365949864064, -0.677336655782164, -0.706075250172249, -0.718165409180243, -0.753647332735972, -0.773208334250904, -0.789018704123706, -0.806061085002831, -0.816875114926754)
# v_upper <- c(-0.521472202180516, -0.444964438492518, -0.327862977696394, -0.235044091838281, 0.406867267767451, 3.21014641683317, 3.48490803223132, 3.82648584443488, 4.07272009850315)
# for (alpha_i in 1:length(alpha_list)){
#     alpha <- alpha_list[alpha_i]
#     a_n <- c(0.2462, 5.020,-0.599)
#     model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=2, 'manual'=FALSE)
#     e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
#     a <- e$Parameter_Limits
#     print(paste(a[1] , v_lower[alpha_i], a[2] , v_upper[alpha_i], sep=" "))
#      
# }
# sink(NULL)
# stop()
control <- list("ncores"=2,'lr' = 0.75,'maxiters' = c(10,10),'halfmax' = 5,'epsilon' = 1e-4,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=0, 'ties'='breslow','double_step'=1, 'guesses'=10)
#v_lower <- c(0.102243677961472, -0.0366388051711791, -0.219114008554502, -0.344171754216377, -1.95449016972488, -2.02344865755178, -2.08312767824572, -2.15183918342603, -2.19804353277016)
#v_upper <- c(0.36980987063755, 0.458228372469313, .539265192368337, 0.575303932769727, 1.07110681914306, 0.756984966349814, 0.815541033177481, 0.88254914662706, .927406355056622)
#for (alpha_i in 1:length(alpha_list)){
#    alpha <- alpha_list[alpha_i]
#    a_n <- c(0.2462, 5.020,-0.599)
#    model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=0, 'manual'=TRUE)
#    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
#    a <- e$Parameter_Limits
#    b <- e$Likelihood_Boundary
#    c <- e$Likelihood_Goal
#    d <- e$Negative_Limit_Found
#    print(paste(a[1], a[2], sep=" "))
#     
#}
# #
start <- Sys.time()
v_lower <- c(4.97252283668956, 4.9349945105648, 4.89804715665926, 4.88084912208962, 4.82369762341988, 4.78721237571926, 4.7546530342797, 4.71603055250556, 4.68938287303871)
v_upper <- c(5.06762896572498, 5.10561529697034, 5.14327069556976, 5.16088604918614, 5.21982880792394, 5.2577860471215, 5.29187760184654, 5.33258757872226, 5.36084782852899)
for (alpha_i in 1:length(alpha_list)){
    alpha <- alpha_list[alpha_i]
    a_n <- c(0.2462, 5.020,-0.599)
    model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=1, 'manual'=TRUE)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
    a <- e$Parameter_Limits
    print(paste(a[1] - v_lower[alpha_i], a[2] - v_upper[alpha_i], sep=" "))
     
}
# #
v_lower <- c(-0.643365949558998, -0.677336655540846, -0.706075250211414, -0.718165409196492, -0.753647332793819, -0.773208334303991, -0.789018704115451, -0.806061085000755, -0.816875114954096)
v_upper <- c(-0.521472203247917, -0.444964438813732, -0.327862977142017, -0.235044092073815, 2.91573713669059, 3.21014641617297, 3.48490803194128, 3.82648584413642, 4.07272009904963)
for (alpha_i in 1:length(alpha_list)){
    alpha <- alpha_list[alpha_i]
    a_n <- c(0.2462, 5.020,-0.599)
    model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=2, 'manual'=TRUE)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
    a <- e$Parameter_Limits
    print(paste(a[1] - v_lower[alpha_i], a[2] - v_upper[alpha_i], sep=" "))
     
}

sink(NULL)
end <- Sys.time()
print(end-start)
stop()