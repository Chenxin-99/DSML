rm(list=ls())
setwd(code_path)
source("CrossValidation.R")
library("caTools")

#for imbalance dataset
n <- 1000
m <- floor(n*400)
score1 <- data.frame(value=rnorm(n, mean=1, sd=1))
score2 <- data.frame(value=rnorm(m, mean=-1, sd=1))
score <- rbind(score1, score2)
label1 <- data.frame(value=rep(1,n))
label2 <- data.frame(value=rep(0,m))
label <- rbind(label1, label2)
data_ROC_n <- data.frame(prob=score$value, obs=label$value)
data_ROC_n <- data_ROC_n[order(-data_ROC_n$prob),]
FTP <- Get_fpr_tpr_pre(data_ROC_n)

#for balance dataset
m2 = floor(n)
score1 = data.frame(value=rnorm(n, mean=1, sd=1))
score2 = data.frame(value=rnorm(m2, mean=-1, sd=1))
score = rbind(score1, score2)
label1 = data.frame(value=rep(1,n))
label2 = data.frame(value=rep(0,m2))
label = rbind(label1, label2)
data_ROC_n = data.frame(prob=score$value, obs=label$value)
data_ROC_n <- data_ROC_n[order(-data_ROC_n$prob),]
FTP2 <- Get_fpr_tpr_pre(data_ROC_n)

#calculate AUC and AUPR
AUC1 <- trapz(FTP$fpr_n, FTP$tpr_n)
AUPR1 <- trapz(FTP$tpr_n, FTP$pre_n)
AUC2 <- trapz(FTP2$fpr_n, FTP2$tpr_n)
AUPR2 <- trapz(FTP2$tpr_n, FTP2$pre_n)

#plot ROC and PR curves
plot(FTP$fpr_n, FTP$tpr_n,type='l',ylab = "TPR",xlab = "FPR",lwd=3)
lines(FTP2$fpr_n, FTP2$tpr_n, col='red')
abline(a=0,b=1,lty=2)
legend("bottomright",inset=0.05,c(paste("imbalance dataset:",round(AUC1,4)),paste("balance dataset:",round(AUC2,4))),col=c("black","red"),lty=1,lwd=3)

plot(FTP$tpr_n, FTP$pre_n,type='l',ylab = "Precision",xlab = "Recall",lwd=3)
lines(FTP2$tpr_n, FTP2$pre_n, col='red')
abline(a=0,b=1,lty=2)
legend("topright",inset=0.05,c(paste("imbalance dataset:",round(AUPR1,4)),paste("balance dataset:",round(AUPR2,4))),col=c("black","red"),lty=1,lwd=3)

