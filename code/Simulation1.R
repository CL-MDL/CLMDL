rm(list=ls())
source("Util_Functions.R")
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(24)
############# Code for simulation 1 in the manuscript


### All simulation settings considered in Simulation 1 
S_sets <- c(6^2, 8^2, 10^2) # Spatial dimension S
TT_sets <- c(100, 200) # Temporal dimension T
delta_sets <- rbind(c(0,0), c(2,0), c(3,0), c(0,6), c(0,10), c(2,2), c(3,3))/10 # Change size


### Select a specific setting in Simulation 1 to run
# i.e. combination of (S, T, delta)
S <- S_sets[1]
TT <- TT_sets[1]
T1 <- T2 <- TT/2 # segment length (change-point is at 0.5*TT)
delta <- delta_sets[3,]

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
final_result <- foreach(rep_index=1:1000, .packages=c('mvtnorm')) %dopar% {
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
  
  pelt_result
}

# Analysis (i.e. check the proportion of experiments where correct number of cp is detected)
sum(sapply(final_result, length)==3)/1000
