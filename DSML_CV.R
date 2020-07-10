rm(list=ls())

setwd(data_path)
RDS <- read.table("DD_MAT.txt") #known effective combination

#drug anatomical therapeutic similarity calculation
SDA1 <- read.table("DD_sim_ATC1.txt")
SDA2 <- read.table("DD_sim_ATC2.txt")
SDA3 <- read.table("DD_sim_ATC3.txt")
SDA1 <- as.matrix(SDA1)
SDA2 <- as.matrix(SDA2)
SDA3 <- as.matrix(SDA3)
SDA <- (SDA1+SDA2+SDA3)/3

RDG <- read.table("DG_MAT.txt") #known drug target
RGG <- read.table("GG.txt")   #known PPI
RDG <- as.matrix(RDG)
RDS <- as.matrix(RDS)
RGG <- as.matrix(RGG)

drug_num <- dim(RDS)[1]
gene_num <- dim(RGG)[1]

library("caTools")

source("Lapla.R")
source("CrossValidation.R")
source("PathSim.R")

#each protein is considered to have self-relation
for (i in 1:dim(RGG)[1]) {
  RGG[i,i] <- 1
}

#calculate protein topoloty similarity
SGG1 <- PathSim(RGG,2)
SGG2 <- PathSim(RGG,3)
SGG3 <- PathSim(RGG,4)
SGG <- (1/3)*(SGG1+SGG2+SGG3)

#runs=1 for leave-one-out cross-validation; runs=10 for five-fold cross-validation
runs <- 1 
drug_num <- dim(RDS)[1]
#record position of all known combinations
p_in_all <- P_positive(RDS)
K <- dim(p_in_all)[1] # if K=5, execute five-fold cross validation; if K=dim(p_in_all)[1], execute leave-one-out cross-validation
num_te <- floor(dim(p_in_all)[1]/K)

TPR_ALL_N <- matrix(nrow=(drug_num*(drug_num-1)/2-(dim(p_in_all)[1]-num_te))+1,ncol=runs*floor(K))
FPR_ALL_N <- matrix(nrow=(drug_num*(drug_num-1)/2-(dim(p_in_all)[1]-num_te))+1,ncol=runs*floor(K))
PRE_ALL_N <- matrix(nrow=(drug_num*(drug_num-1)/2-(dim(p_in_all)[1]-num_te))+1,ncol=runs*floor(K))

#set parameter combination
lamdaS <- 0.1
lamdaA <- 0.1
lamdaT <- 0.7


for (r in 1:runs) {
  
  RDS_v <- RDS
  
  tep_pos_set <- sample(dim(p_in_all)[1],dim(p_in_all)[1])
  tep_pos_set_all <- rep(0,num_te)
  
  for (i in 1:floor(K)) {
    t_p <- 1
    RDS_v <- RDS
    for (j in ((i-1)*num_te+1):(i*num_te)) {
      RDS_v[p_in_all[tep_pos_set[j],1],p_in_all[tep_pos_set[j],2]] <- 0
      RDS_v[p_in_all[tep_pos_set[j],2],p_in_all[tep_pos_set[j],1]] <- 0
    }
    
    #obtain Laplacian matrix
    LDG <- Lapla_nor(SGG)
    LDS <- Lapla_nor(RDS_v)
    LDA <- Lapla_nor(SDA)

    #obtain normalized Laplacian matrix
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

    rate <- 0.01
    X2 <- X1-rate*gradient_x
    Y2 <- Y1-rate*gradient_y
    iter <- 1

    while (sqrt(sum(X2-X1)^2)>0.01 & iter<100) {
      X1 <- X2
      Y1 <- Y2
      gradient_x <- 2*((1+gama)*X1+lamdaS*LLDS%*%X1+lamdaA*LLDA%*%X1-gama*t(Y1)%*%Y1-X_R)
      gradient_y <- 2*(Y1+lamdaT*LLDT%*%Y1-2*gama*Y1%*%X1+2*gama*Y1%*%t(Y1)%*%Y1-Y_R)

      X2 <- X1-rate*gradient_x
      Y2 <- Y1-rate*gradient_y

      print(sqrt(sum(X2-X1)^2))
      iter <- iter + 1
    }
    F <- X2+t(X2)
    
    data_ROC_n <- Get_Test_Score(F,RDS_v,RDS)
    FTP <- Get_fpr_tpr_pre(data_ROC_n)
    TPR_ALL_N[,((r-1)*floor(K)+i)] <- FTP$tpr_n
    FPR_ALL_N[,((r-1)*floor(K)+i)] <- FTP$fpr_n
    PRE_ALL_N[,((r-1)*floor(K)+i)] <- FTP$pre_n
  }
}

tpr_p <- rowMeans(TPR_ALL_N)
fpr_p <- rowMeans(FPR_ALL_N)
pre_p <- rowMeans(PRE_ALL_N)

AUC <- trapz(fpr_p,tpr_p)
print(AUC)














