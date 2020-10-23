library(swCRTdesign)
library(varhandle)
library(lme4)
library(lmerTest)
library("argparse")
source("Imb.r")
source("paperfuns.R")

## -----------------------------------------
## load command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--I", type = "double",default = 24,
                    help = "number of clusters")
parser$add_argument("--tau", type = "double", default = 0.3,
                    help = "random cluster effect")
parser$add_argument("--eta", type = "double", default = 0.1,
                    help = "random treatment effect")
parser$add_argument("--rho", type = "double", default = 0.1,
                    help = "random cluster-treatment correlatioin")
parser$add_argument("--theta", type = "double", default = 0.15,
                    help = "treatment effect under the alternative")
parser$add_argument("--nsim", type = "double",default = 1000,
                    help = "number of simulations per setting")
parser$add_argument("--nsim-perjob", type = "double",default = 10,
                    help = "number of simulations per job")
args <- parser$parse_args()



## -----------------------------------------
## set up the values of necessary parameters
## -----------------------------------------
I = args$I
J = 5 
N = 100 #no. of subjects per cluster-time unit
tau = args$tau
eta = args$eta
rho = args$rho
gamma = 0.1
swDesign = swDsn( clusters=rep(I/(J-1),J-1))
sigma = 1
#dummy representation has one less column than the number of factors
#so omega has 1+2+1=4 entries
omega = matrix(c(1,1,2,1),nrow=4,ncol=1) #coefficient of cluster-level covariates
Theta = c(0,args$theta)
opt=2 #type of imbalance score to be used 1:trt-ctl balance 2: mean balance


## -----------------------------------------
## set up a grid of parameters to cycle over
## -----------------------------------------
njob_percombo = args$nsim / args$nsim_perjob
param_grid <- expand.grid(sim_id = 1:njob_percombo,Theta=Theta)

## -----------------------------------------
## get current dynamic arguments
## -----------------------------------------
## get job id from scheduler
job_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

## current args
current_dynamic_args <- param_grid[job_id, ]
theta = current_dynamic_args$Theta

Seq <- swDesign$swDsn
#create candidate set (delete replicates from 20000 permutations)
cand <- Cand_list(I,J)


## -----------------------------------------
## run the simulation  nsim_perjob times 
## -----------------------------------------
# set current seed
current_seed <- job_id*123 + 12345 #default seed
set.seed(current_seed)
result <- matrix(nrow = 1, ncol = 3)

for (t in 1:args$nsim_perjob) {
  #one continuous covariate, one categorical factor with level 2, 
  #one categorical factor with level 3.
  Z1 <- matrix(data = runif(I,min=-1,max=1), nrow = I, ncol = 1)
  Z2 <- matrix(data = t(rmultinom(I, size = 1, prob = c(1/3,1/3,1/3))),
               nrow = I, ncol = 3)
  Z3 <- matrix(data = t(rmultinom(I, size = 1, prob = c(0.5,0.5))),
               nrow = I, ncol = 2)
  Z_init = cbind(Z1,Z2,Z3)
  #get the weight vector
  W = c(1,apply(Z2==1,2,sum)/I,apply(Z3==1,2,sum)/I)
  
  cand_final <- Final_list(cand=cand,Z_init=Z_init,
                           p=1,
                           opt=opt,Seq=Seq,W=W)
  n_seq <- dim(cand_final)[2] #number of sequences in the final candidate set
  #drop the first column of each categorical factors to avoid rank deficiency
  Z_init <- Z_init[,-c(2,5)]
  
  #Generate the outcome Y, select the treatment allocation for 
  #level of constraint q varying from 0.1 to 1
  for (q in seq(0.1,1,0.1)) {
    #Generate outcome data
    Y <- swSim(
      swDesign,
      family = "gaussian",
      n = N,
      mu0 = 0,
      mu1 = theta,
      time.effect = (0:(J - 1)),
      sigma = sigma,
      tau = tau,
      eta = eta,
      rho = rho,
      gamma = gamma
    )
    #randomly select a sequence from final candidate set with idx
    # constrained within the first q-th percentile
    idx <- sample(1:floor(n_seq*q))[1]
    candi <- cand_final[, idx]
    seq_Z <- mat.sort(as.matrix(cbind(candi, Z_init)), "candi")
    Z <-  as.matrix(seq_Z[,-which(colnames(seq_Z)
                                  %in% c("candi"))])
    Z_all <- apply(Z, 2, rep, each = J * N)
    colnames(Z_all) <- (1:ncol(Z_all))
    Y$Z_all <- Z_all
    Y$response.final <- as.vector(Y$response.var + Z_all %*% omega)
    
    Y_sw <- swSummary(response.final,tx.var,time.var,cluster.var,
                      Y,type="mean")
    Ymean <- Y_sw$response.cluster
    swmat <- swDesign$swDsn
    
    #fit model
    if (eta == 0) {
      fit_unadj <- lmer(
        response.final ~ time.var + tx.var
        + (1 |
             cluster.var) 
        + (1 | cluster.var:time.var)  ,
        data = Y
      )
      fit_adj <- lmer(
        response.final ~ time.var + tx.var + Z_all
        + (1 |
             cluster.var) 
        + (1 | cluster.var:time.var)  ,
        data = Y
      )
    } else{
      fit_unadj <- lmer(
        response.final ~ time.var + tx.var +
          (1 | cluster.var:time.var)  + (tx.var |
                                           cluster.var),
        data = Y
      )
      fit_adj <- lmer(
        response.final ~ time.var + tx.var + Z_all +
          (1 | cluster.var:time.var)  + (tx.var |
                                           cluster.var),
        data = Y
      )
    }
    
    ###########################################
    #robust permutstion inference(unstratified)
    ###########################################
    thetahat_unadj <- ubest(Ymean,swmat)
    Z_unadj <- thetahat_unadj/sqrt(vardpc(Ymean,swmat,0))
    pval_unadj <- pnorm(-abs(Z_unadj)) * 2
    
    result <- rbind(result,summary(fit_unadj)$coefficients["tx.var",
                                                           c("Estimate","Std. Error","Pr(>|t|)") ])
    result <- rbind(result,summary(fit_adj)$coefficients["tx.var",
                                                         c("Estimate","Std. Error","Pr(>|t|)") ])
    result <- rbind(result,c(thetahat_unadj,sqrt(vardpc(Ymean,swmat,0)),
                             pval_unadj))
    
    cat("current iteration is", t,",with level of constraint", q, "\n")
  }
  
}

## ------------------------------------------------------
## Organize the result, save them as Rdata for later use 
## ------------------------------------------------------
num_p = 10
result = result[-1,]
result <- as.data.frame(result)
result$Reject = as.numeric(result$`Pr(>|t|)` <= 0.05)
result$Method = rep(c("lme-Unadjusted","lme-Adjusted","perm-unstratified"),
                    args$nsim_perjob*num_p)
result$p = rep(rep(seq(0.1,1,0.1),each=3),args$nsim_perjob)
result$I = I
result$tau = tau
result$eta = eta
result$theta = theta
file_name <- paste0("output/","I=", I,"tau=",tau, "eta=",eta,
                    "theta=",theta,"bscore=",opt, "job_no=", 
                    current_dynamic_args$sim_id,".rds")
saveRDS(result, file = file_name)


