rm(list=ls())
source("Util_Functions.R")
library(mvtnorm)
set.seed(1234)

############# DGP for simulating a spatio-temporal process ############
### Number of CP
cp_sets <- c(0,1,2)
for(cp in cp_sets){
  ### Spatial dimension S
  S <- 6^2
  min.dist <- 1
  max.lati <- max.long <- sqrt(S)
  
  n.lati <- max.lati/min.dist
  n.long <- max.long/min.dist
  
  is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x - round(x)) < tol # function to check if integer within small tol
  if(!is.wholenumber(n.lati) || !is.wholenumber(n.long)) stop("Not divisible n.lati or n.long!!!")
  
  ### Generate matrix of 2D regular grid and distance matrix
  mat.lattice <- cbind(rep(seq(1, max.lati, by=min.dist), n.long), rep(seq(1, max.long, by=min.dist), each=n.lati))
  S.dist <- round(as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE)), 10) # round to 10 decimal place to remove extra distance combinations
  
  ### Time dimension T and CP location
  # No change-point
  if(cp==0){
    data_name <- 'y0'
    TT <- 100
    theta1 <- c(-0.5, 0.6, 1, 0) # c(phi, rho, sigma2, mu)
    y <- sim.y(theta=theta1, S.dist=S.dist, TT=TT, T.burn=100)
  }
  # One change-point
  if(cp==1){
    data_name <- 'y1'
    T1 <- T2 <- 50 # segment length
    TT <- T1 + T2 # total time
    
    # True model parameters and change size
    delta <- 0.3
    theta1 <- c(-0.5, 0.6, 1, 0) # c(phi, rho, sigma2, mu)
    theta2 <- c(-0.5 + delta, 0.6 + delta, 1, 0)
    
    # Simulate a dataset with one single change-point at 0.5*T
    y <- rbind(sim.y(theta=theta1, S.dist=S.dist, TT=T1, T.burn=100),
               sim.y(theta=theta2, S.dist=S.dist, TT=T2, T.burn=100))
  }
  # Two change-point
  if(cp==2){
    data_name <- 'y2'
    # Time dimension T and CP location
    T1 <- T2 <- 33
    T3 <- 34
    TT <- T1 + T2 + T3
    
    # True model parameters and change size
    delta <- 0.2
    theta1 <- c(-0.5, 0.6, 1, 0) # c(phi, rho, sigma2, mu)
    theta2 <- c(-0.5 + delta, 0.6 + delta, 1, delta)
    
    # Simulate data
    y <- rbind(sim.y(theta=theta1, S.dist, T=T1, T.burn=100),
               sim.y(theta=theta2, S.dist, T=T2, T.burn=100),
               sim.y(theta=theta1, S.dist, T=T3, T.burn=100))
  }
  
  data <- list(y=y, S.dist=S.dist)
  
  saveRDS(data, file=paste0(data_name, '.RDS'))
}

