rm(list=ls())

setwd(data_path)
RDS <- read.table("DD_MAT.txt")
SDA1 <- read.table("DD_sim_ATC1.txt")
SDA2 <- read.table("DD_sim_ATC2.txt")
SDA3 <- read.table("DD_sim_ATC3.txt")
RDG <- read.table("DG_MAT.txt")
RGG <- read.table("GG.txt")
drug_name <- read.table("drug_ALL.txt")

RDG <- as.matrix(RDG)
RDS <- as.matrix(RDS)
RGG <- as.matrix(RGG)
SDA1 <- as.matrix(SDA1)
SDA2 <- as.matrix(SDA2)
SDA3 <- as.matrix(SDA3)
SDA <- (SDA1+SDA2+SDA3)/3

#each protein is considered to have self-relation
for (i in 1:dim(RGG)[1]) {
  RGG[i,i] <- 1
}

drug_num <- dim(RDS)[1]
gene_num <- dim(RGG)[1]

source("BMA.R")
source("Lapla.R")
source("CrossValidation.R")
source("PathSim.R")

for (i in 1:dim(RGG)[1]) {
  RGG[i,i] <- 0
}

SGG1 <- PathSim(RGG,2)
SGG2 <- PathSim(RGG,3)
SGG3 <- PathSim(RGG,4)
SGG <- (1/3)*(SGG1+SGG2+SGG3)

lamdaA <- 0.1
lamdaS <- 0.1
lamdaT <- 0.7

RDS_v <- RDS

LDG <- Lapla_nor(SGG)  
LDS <- Lapla_nor(RDS_v)
LDA <- Lapla_nor(SDA)
LDC <- Lapla_nor(SDC)

LLDT <- Lapla_matrix_I(LDG)
LLDS <- Lapla_matrix_I(LDS)
LLDA <- Lapla_matrix_I(LDA)

b <- 1
I_D <- diag(drug_num)
I_G <- diag(gene_num)

X_R <- RDS_v
Y_R <- t(RDG)
X1 <- X_R
Y1 <- Y_R

gama <- 0.001
#calculate gradient
gradient_x <- 2*((1+gama)*X1+lamdaS*LLDS%*%X1+lamdaA*LLDA%*%X1-gama*t(Y1)%*%Y1-X_R)
gradient_y <- 2*(Y1+lamdaT*LLDT%*%Y1-2*gama*Y1%*%X1+2*gama*Y1%*%t(Y1)%*%Y1-Y_R)

rate <- 0.001
X2 <- X1-rate*gradient_x
Y2 <- Y1-rate*gradient_y

iter <- 1
while (sqrt(sum(X2-X1)^2)>0.01) {
  X1 <- X2
  Y1 <- Y2
  #calculate gradient
  gradient_x <- 2*((1+gama)*X1+lamdaS*LLDS%*%X1+lamdaA*LLDA%*%X1-gama*t(Y1)%*%Y1-X_R)
  gradient_y <- 2*(Y1+lamdaT*LLDT%*%Y1-2*gama*Y1%*%X1+2*gama*Y1%*%t(Y1)%*%Y1-Y_R)

  X2 <- X1-rate*gradient_x
  Y2 <- Y1-rate*gradient_y
  print(sqrt(sum(X2-X1)^2))
  iter <- iter + 1
}
F <- X2+t(X2)

#rank potential drug pairs with synergy score
data_ROC_n <- Get_Test_Score_List(F,RDS_v,RDS)
drug_com <- data.frame(matrix(0,nrow = 100,ncol = 2))
for (i in 1:100) {
  data_out <- data.frame(drug_name[data_ROC_n$r[i],1], drug_name[data_ROC_n$c[i],1])
  write.table(data_out,file="ALL_Rank3.txt",append =TRUE)
}
