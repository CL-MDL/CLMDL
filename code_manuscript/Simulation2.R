rm(list=ls())
source("Util_Functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(24)
############# Code for simulation 2 in the manuscript


### All simulation settings considered in Simulation 2
rep_times <- 100
S_sets <- c(6^2, 8^2, 10^2) # Spatial dimension S
TT_sets <- c(200) # Temporal dimension T
delta_sets <- rbind(c(2,0), c(3,0), c(2,2))/10 # Change size


### Select a specific setting in Simulation 2
# i.e. combination of (S, T, delta)
S <- S_sets[1]
TT <- TT_sets[1]
T1 <- T2 <- TT/2 # segment length (change-point is at 0.5*TT)
delta <- delta_sets[2,]

# True model parameters in the first and second segment c(phi, rho, sigma2)
theta1 <- c(-0.5, 0.6, 1)
theta2 <- c(-0.5 + delta[1], 0.6 + delta[2], 1)


### Generate matrix of 2D regular grid and distance matrix
min.dist <- 1
max.lati <- max.long <- sqrt(S)
n.lati <- max.lati/min.dist
n.long <- max.long/min.dist
mat.lattice <- cbind(rep(seq(1, max.lati, by=min.dist), n.long), rep(seq(1, max.long, by=min.dist), each=n.lati))
S.dist <- round(as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE)), 10) # round to 10 decimal place to remove extra distance combinations


### Parameters for implementing CLDML
# Set spatial and temporal lag in computing CL
s.lag <- t.lag <- 1

# Collect all pairs of observations within s.lag distance (will be used in PL)
s.comb <- s.dist.pair(S.dist, s.lag) # time.lag > 0
s.comb.0 <- s.dist.pair.0(S.dist, s.lag) # time.lag = 0
s.len <- unique(s.comb[,3])
s.len <- s.len[s.len>0] # all unique spatial distance


######### Main function ##########
final_result <- foreach(rep_index=1:rep_times, .combine='rbind', .packages='mvtnorm') %dopar% {
  # Simulate a dataset with one single change-point at 0.5*T
  y <- rbind(sim.y(theta=theta1, S.dist=S.dist, TT=T1, T.burn=100),
             sim.y(theta=theta2, S.dist=S.dist, TT=T2, T.burn=100))
  
  # Calculate average number of times an obs used in CL for Ck
  pair_stat <- D.cal(y, S.dist, s.lag, t.lag)
  Ck <- mean(pair_stat$use.obs)
  remedy <- pair_stat$remedy # number of marginal likelihood needed for correcting edge effect
  
  # Calculate K for pruning step
  p.min <- p.max <- length(theta1)
  K <- floor(Ck)*(floor(log(S*TT))*(p.min/2 - p.max) + (2 + p.max)*log(2) - floor(log(TT)))
  
  # initial parameter estimation based on entire data using the 4-parameter spatial auto-regressive model
  ini <- optim(c(0,1,1), pl, y=y, S.dist=S.dist, t.lag=t.lag, remedy=remedy, mu_zero=T,
               method="L-BFGS-B", lower=c(-0.7,0.1,1e-3), upper=c(0.7,3,5))$par
  
  # Run PELT for minimizing CLMDL (this is the default algorithm in the paper)
  pelt_result <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K,
                      t.lag=t.lag, min.length=0.2*TT, ini=ini, res=2, mu_zero=T)
  
  # Construct CI for the estimated change-points
  num_cp <- length(pelt_result)-2
  if(num_cp==1){
    est_cp <- pelt_result[2]
    # Parameter estimation on the pre-change segment
    y.seg1 <- y[1:est_cp,]
    pmle1 <- optim(ini, pl, y=y.seg1, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
                   method="L-BFGS-B", lower=c(-0.7,0.1,1e-3), upper=c(0.7,3,5), ts=F, mu_zero=T)
    # Parameter estimation on the post-change segment
    y.seg2 <- y[(est_cp+1):TT,]
    pmle2 <- optim(ini, pl, y=y.seg2, S.dist=S.dist, t.lag=t.lag, remedy=remedy,
                   method="L-BFGS-B", lower=c(-0.7,0.1,1e-3), upper=c(0.7,3,5), ts=F, mu_zero=T)
    
    seg1_tmp <- nrow(y.seg1)
    seg2_tmp <- nrow(y.seg2)
    TT_tmp <- seg1_tmp + seg2_tmp
    
    # CI construction
    B <- 100
    theta1_est <- pmle1$par
    theta2_est <- pmle2$par
    q_emp <- c()
    for(b in 1:B){
      yB <- rbind(sim.y(theta=theta1_est, S.dist=S.dist, TT=seg1_tmp, T.burn=100),
                  sim.y(theta=theta2_est, S.dist=S.dist, TT=seg2_tmp, T.burn=100))
      
      q_emp[b] <- Q.grid(lambda.hat=seg1_tmp, TT=TT_tmp, S.dist=S.dist, y.full=yB,
                         theta1=theta1_est, theta2=theta2_est,
                         t.lag=t.lag, remedy=remedy, q.bound=round(0.2*TT_tmp), ts=F, mu_zero=T)
    }

    ci_lower <- -ceiling(-quantile(q_emp, probs=c(0.05,0.025,0.005)))
    ci_upper <- ceiling(quantile(q_emp, probs=c(0.95,0.975,0.995)))

    ci_result90 <- est_cp+c(ci_lower[1], ci_upper[1])+c(-1,1)
    ci_result95 <- est_cp+c(ci_lower[2], ci_upper[2])+c(-1,1)
    ci_result99 <- est_cp+c(ci_lower[3], ci_upper[3])+c(-1,1)
  }
  
  sum_result <- list(est=pelt_result, ci90=ci_result90, ci95=ci_result95, ci99=ci_result99, q_emp=q_emp)
  sum_result
}

save.image("Simulation2.RData")

### Analysis (i.e. check change-point estimation accuracy and check CI coverage)
num_cp <- cp <- coverage <- c()
for(b in 1:rep_times){
  tmp90 <- final_result[b,]$ci90-TT/2
  tmp95 <- final_result[b,]$ci95-TT/2
  tmp99 <- final_result[b,]$ci99-TT/2
  num_cp <- c(num_cp, length(final_result[b,]$est))
  cp <- rbind(cp, final_result[b,]$est)
  
  tmp_cover <- c(tmp90[1]<=0 & tmp90[2]>=0, tmp95[1]<=0 & tmp95[2]>=0, tmp99[1]<=0 & tmp99[2]>=0)
  coverage <- rbind(coverage, tmp_cover)
}

# estimation accuracy
sum(num_cp==3)
apply(cp, 2, mean)/TT
apply(cp, 2, sd)/TT

# CI coverage
apply(coverage, 2, mean)

