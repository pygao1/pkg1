########################
# Functions
########################
library(tidyverse)
require(data.table)

##############################################
# Compute the imbalance score(modification of the B score proposed by 
#Raab and Butcher for continuous covariate) 
Imb_cont1 = function(Z,Seq){
  # Z: matrix of covariates, each covariate is organized as seperate column
  # Seq: SW sequence of control(0)-intervention(1) for each cluster
  Z_std <- scale(Z)
  ti <- apply(Seq==1,1,which.max)
  sum(apply((2*ti -ncol(Seq) -2) * Z_std,2,sum)^2)
}

# treatment-control balance score, with pre-assigning weight of each
# covariate (1 for continuous, relative frequency for each level of
# categorical covariates)
Imb_tc = function(Z,Seq,W){
  Z_std <- scale(Z)
  ti <- apply(Seq==1,1,which.max)
  score = ((2*ti -ncol(Seq) -2) %*% Z)^2 %*% W
  return(score)
}

# Mean balance score for continuous covariate
Imb_cont2 = function(Z,Seq){
  # Z: matrix of covariates, each covariate is organized as seperate column
  # Seq: SW sequence of control(0)-intervention(1) for each cluster
  #standardize the covariate values
  Z_std <- scale(Z)
  # find the time wave index of each cluster
  Ti = apply(Seq==1,1,which.max)
  Ti = as.factor(Ti)
  n_waves = length(levels(Ti))
  # Get the mean values under eachh time waves
  Z_iw = matrix(0,nrow=n_waves,ncol=ncol(Z_init))
  for (w in 1:n_waves){
    Z_iw[w,] = apply(matrix(Z_std[Ti==levels(Ti)[w],],ncol=ncol(Z_init)),
                     2,mean)
  }
  return(sum(Z_iw^2))
}

# mean balance score, with pre-assigning weight of each
# covariate (1 for continuous, relative frequency for each level of
# categorical covariates)
Imb_mean = function(Z,Seq,W){
  # Z: matrix of covariates, each covariate is organized as seperate column
  # Seq: SW sequence of control(0)-intervention(1) for each cluster
  #standardize the covariate values
  Z_std <- scale(Z)
  # find the time wave index of each cluster
  Ti = apply(Seq==1,1,which.max)
  Ti = as.factor(Ti)
  n_waves = length(levels(Ti))
  # Get the mean values under eachh time waves
  Z_iw = matrix(0,nrow=n_waves,ncol=ncol(Z_init))
  for (w in 1:n_waves){
    Z_iw[w,] = apply(matrix(Z_std[Ti==levels(Ti)[w],],ncol=ncol(Z_init)),
                     2,mean)
  }
  return(apply(Z_iw^2,2,sum) %*% W)
}


Imb_mean2 = function(Z,Seq,W){
  # Z: matrix of covariates, each covariate is organized as seperate column
  # Seq: SW sequence of control(0)-intervention(1) for each cluster
  #standardize the covariate values
  Z_std <- scale(Z)
  # find the time wave index of each cluster
  Ti = as.factor(apply(Seq==1,1,which.max))
  n_c = nrow(Z)
  n_waves = length(levels(Ti))
  n_rep = n_c / n_waves
  # Get the mean values under eachh time waves
  Z_iw = matrix(0,nrow=n_waves,ncol=ncol(Z_init))
  for (w in 1:n_waves){
    Z_iw[w,] = apply(Z_std[(n_rep*(w-1)+1):(n_rep*w),],
                     2,mean)
  }
  return(apply(Z_iw^2,2,sum) %*% W)
}

##############################################
Cand_list = function(I,J,cand_size=20000){
  # create candidate set(with size=cand_size) of possible intervention sequences, 
  # represented as the last time before switching to treatment
  # dubplicate sequence is deleted
  # I/(J-1) is assumed to be an interger so there're integer number of complete sequences
  cand <- matrix(0,nrow=cand_size,ncol=I)
  cand_null <- rep(1:(J-1),each=I/(J-1))
  for (t in 1:cand_size){
    cand[t,] <- sample(cand_null,I)
  }
  #cand <- as.data.frame(cand)
  #cand <- distinct(cand)
  cand <- t(unique(cand))
  #cand <- as.matrix(t(cand))
  cand
}

##############################################
#sort a matrix according to a specified column
mat.sort <- function(mat,sort,decreasing=FALSE){
  m<-do.call("order",c(as.data.frame(mat[,sort]),decreasing=decreasing))
  mat[m,]
} 



#################################
Final_list = function(cand,Z_init,p,opt=1,Seq,W){
  # the function to select the final candidate set 
  # with smallest imbalance score
  # cand: candidate list returned by the function Cand_list
  # Z_init:covariate vector(without permutation) with size 
  # I(no.of cluster)*S(no.of covariate)
  # p: proportion of candidate set to be selected
  # opt: if 1, use Imb_cont1, if 2, use Imb_cont2
  # Seq: SW sequence of control(0)-intervention(1) for each cluster
  #W : weight of each covqriate and levels of categorical factors
  size0 <- dim(cand)[2] #initial size of candidate set
  ImbScore <- rep(NA,size0)
  for (k in 1:size0){
    candk <- as.vector(cand[,k])
    seq_Z <-as.data.frame(mat.sort(cbind(candk,Z_init),"candk"))
    Z <- as.matrix(seq_Z[ , -which(names(seq_Z) 
    %in% c("candk"))])
    if (opt==1){
    ImbScore[k] = Imb_tc(Z,Seq,W)} else if (opt==2){
      ImbScore[k] = Imb_mean2(Z,Seq,W)
    }
  }
  score_list <- mat.sort(as.matrix(cbind(score=ImbScore,index=1:size0)),"score")
  size1 <- floor(size0*p) #final set of candidate set
  final_index <- score_list[1:size1,"index"]
  cand_final <- cand[,final_index]
  cand_final
}
