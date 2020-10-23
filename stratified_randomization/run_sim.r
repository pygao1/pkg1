###################################################################
#Simulation for stratified&simple randomization with model and 
#permutation-based inference
###################################################################
library(swCRTdesign)
library(varhandle)
library(lme4)
library(lmerTest)
library("argparse")
source("paperfuns.R")

## -----------------------------------------
## load any command line arguments
parser <- ArgumentParser()
parser$add_argument("--I", type = "double",default = 8,
                    help = "number of clusters")
parser$add_argument("--tau", type = "double", default = 0.1,
                    help = "random cluster effect")
parser$add_argument("--eta", type = "double", default = 0,
                    help = "random treatment effect")
parser$add_argument("--rho", type = "double", default = 0,
                    help = "random cluster-treatment correlatioin")
parser$add_argument("--theta", type = "double", default = 0.2,
                    help = "treatment effect under the alternative")
parser$add_argument("--nsim", type = "double",default = 5000,
                    help = "number of simulations per setting")
parser$add_argument("--nsim-perjob", type = "double",default = 500,
                    help = "number of simulations per job")
args <- parser$parse_args()
## -----------------------------------------

## #####################################################
## set up the values of other parameters 
## #####################################################
I = args$I
J = 5 
N = 100 #no. of subjects per cluster-time unit
tau = args$tau
eta = args$eta
rho = args$rho
gamma = 0.1
swDesign = swDsn( clusters=rep(I/(J-1),J-1))
sigma = 1
omega = 1
L = 2 #no. of levels of the cluster-level factor
Rand = c('Simple','Stratified')
Theta = c(0,args$theta)


## -----------------------------------------
## set up a grid of parameters to cycle over
## -----------------------------------------
njob_percombo = args$nsim / args$nsim_perjob
param_grid <- expand.grid(sim_id = 1:njob_percombo,Rand = Rand,Theta = Theta)

## -----------------------------------------
## get current dynamic arguments
## -----------------------------------------
## get job id from scheduler
job_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
## current args
current_dynamic_args <- param_grid[job_id, ]
theta = current_dynamic_args$Theta
rand = current_dynamic_args$Rand

## --------------------------------------------------------
## set up the simulation and store thhe result into a table
## --------------------------------------------------------
result <- matrix(nrow = 1, ncol = 3)
set.seed(1000) #default seed
for (t in 1:args$nsim_perjob){
  #Generate outcome data
  Y <- swSim(swDesign,family = "gaussian",n=N,mu0=0,mu1=theta,tau=tau,
             gamma=gamma,eta=eta,rho=rho, 
             time.effect=(0:(J-1)),sigma=sigma)
  
  if (rand == 'Simple'){
    #simple randomization
    #permute covariate value size randomly among all sequences
    size <- sample(rep(0:(L-1),I/L))
    size_all <- rep(size,each=J*N)
  } else if (rand == 'Stratified'){
    #Constrained randomization
    #permute covariate value size within each strata
    size <- vector()
    for (i in 1:(J-1)){ 
      size <- c(size,sample(rep(0:(L-1),I/(J-1)/L)))
    }
    size_all <- rep(size,each=J*N)
  }
  
  Y$response.final <- Y$response.var + size_all * omega
  Y$size_all <- size_all
  
  Y_sw <- swSummary(response.final,tx.var,time.var,cluster.var,
                    Y,type="mean")
  
  Ymean <- Y_sw$response.cluster
  
  swmat <- swDesign$swDsn
  
  ###############
  #model fitting
  ###############
  if (eta==0){
    fit_unadj <- lmer(response.final ~  tx.var 
                    + time.var 
                    +(1|cluster.var) 
                    +(1|cluster.var:time.var)  
                    ,data=Y) 
    fit_adj <- lmer(response.final ~ tx.var + size_all
                  +(1|cluster.var) 
                  +time.var
                  +(1|cluster.var:time.var)  
                  ,data=Y)  
  } else
  {
    fit_unadj <- lmer(response.final ~  tx.var 
                      + time.var 
                      + (1|cluster.var:time.var)  
                      +(tx.var|cluster.var)
                      ,data=Y)
    fit_adj <- lmer(response.final ~  tx.var + size_all 
                    +   (1|cluster.var:time.var)  
                    + time.var 
                    +(tx.var|cluster.var)
                    ,data=Y) 
    }
  #############################
  #robust permutstion inference
  #############################
  #unadjusted
  thetahat_unadj <- ubest(Ymean,swmat)
  Z_unadj <- thetahat_unadj/sqrt(vardpc(Ymean,swmat,0))
  pval_unadj <- pnorm(-abs(Z_unadj)) * 2
  #adjusted
  thetahat_adj <- subest(Ymean,swmat,size+1)
  Z_adj <- thetahat_adj/sqrt(svardpc(Ymean,swmat,0,size+1))
  pval_adj <- pnorm(-abs(Z_adj)) * 2
  
  #add record to the combined dataframe
  result <- rbind(result,summary(fit_unadj)$coefficients["tx.var",
                      c("Estimate","Std. Error","Pr(>|t|)") ])
  result <- rbind(result,summary(fit_adj)$coefficients["tx.var",
                      c("Estimate","Std. Error","Pr(>|t|)") ])
  result <- rbind(result,c(thetahat_unadj,sqrt(vardpc(Ymean,swmat,0)),
                           pval_unadj))
  result <- rbind(result,c(thetahat_adj,sqrt(svardpc(Ymean,swmat,0,size+1)),
                           pval_adj))
  cat("current iteration is", t,"\n")
}

## ------------------------------------------------------
## Organize the result, save them as Rdata for later use 
## ------------------------------------------------------
result = result[-1,]
result <- as.data.frame(result)
result$Rand = rep(rand,args$nsim_perjob*4)
result$I = rep(I,args$nsim_perjob*4)
result$tau = rep(tau,args$nsim_perjob*4)
result$theta = rep(theta,args$nsim_perjob*4)
result$eta = rep(eta,args$nsim_perjob*4)
#result$omega = rep(omega,args$nsim_perjob*4)
result$Reject = as.numeric(result$`Pr(>|t|)` <= 0.05)
result$Method = rep(rep(c("model","perm"),each=2),args$nsim_perjob)
result$covar_adj = rep(c("Unadj","Adj"),args$nsim_perjob*2)
file_name <- paste0("output/","I=", I,"tau=",tau,"eta=",eta,
                    "theta=",theta,"Rand=",rand, "job_no=", 
                    current_dynamic_args$sim_id,".rds")
saveRDS(result, file = file_name)




