library(survival)
library(parallel)
library(dplyr)
library(data.table)
library(Rcpp)
library(RcppEigen)


Rcpp::sourceCpp("Generate.cpp")



#string fname, NumericVector Pb, int nm_xtra_rows, NumericVector Pcs, StringVector RowNames


Pb <- c(70,1)

nm_xtra_rows <- 5
Pc <- c(.1/2000,.015,.1/100,.0001,.0001,.0001)

row_names <- c("ID","YOB","age_entry","sexm","edu","duration","cumulative_dose_lung_lag10")
for (i in 1:nm_xtra_rows){
    row_names <- c(row_names,paste("lung",i,sep=""))
    Pc <- c(Pc, .1*(i+.1)/(nm_xtra_rows))
}
row_names <- c(row_names,"base_hazard","HazardRatio","Lung_Status","Base_Lung_Status")
fname <- "temp.csv"
RunGenerator(fname, Pb, nm_xtra_rows,Pc,row_names,1e6)

#colTypes=c("integer","double","double","double","integer","character","integer","integer","character","integer","integer", "integer","integer","character","character","character","numeric","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer")
df <- fread(fname,nThread=detectCores(),data.table=TRUE,header=TRUE,verbose=TRUE,fill=TRUE)
setkeyv(df, "age_entry")
num_p <- nrow(df)
df=df[age_entry>0,]
df=df[age_entry<150,]
df=df[HazardRatio<100,]
num_a <- nrow(df)

print(num_p-num_a)
print(num_a)

x <- df$age_entry
y <- df$base_hazard

print("plotting risk")
jpeg(paste("risk_plot.jpg",sep=""))
plot(x,y, type="p")
dev.off()

jpeg("hist.png")
hist(df$HazardRatio,10)
dev.off()

dfend <- df[Lung_Status==1, ]
l_num = unlist(unique(df[,ID]),use.names=FALSE)
#
t<- c(0)
n <- c(1)
tu <- unlist(unique(dfend[,age_entry]))
for (i in tu[1]:tu[length(tu)]){
    #
    df0 <- df[age_entry<i,]
    num = length(unlist(unique(df0[,ID]),use.names=FALSE))
    ev = sum(df0[, Lung_Status])
    #
    if (num>0){
        t <- c(t,i)
        temp <- (num - ev)/num
        n <- c(n, temp)
    }
}

jpeg(paste("Surv.jpg",sep=""), units="in", width=5, height=5, res=1200)
plot(t,n,col='red',type='l')
dev.off()
