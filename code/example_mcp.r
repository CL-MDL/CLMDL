rm(list=ls())
source("aux_functions.r")
library(mvtnorm)
set.seed(1234)

### DGP for the spatio-temporal process
# Spatial dimension S
S <- 6^2
min.dist <- 1
max.lati <- max.long <- sqrt(S)

n.lati <- max.lati/min.dist
n.long <- max.long/min.dist

is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x - round(x)) < tol # function to check if integer within small tol
if(!is.wholenumber(n.lati) || !is.wholenumber(n.long)) stop("Not divisible n.lati or n.long!!!")

# Generate matrix of 2D regular grid and distance matrix
mat.lattice <- cbind(rep(seq(1, max.lati, by=min.dist), n.long), rep(seq(1, max.long, by=min.dist), each=n.lati))
S.dist <- round(as.matrix(dist(mat.lattice, upper=TRUE, diag=TRUE)), 10) # round to 10 decimal place to remove extra distance combinations

# Time dimension T and CP location
T1 <- T2 <- T3 <- 50
T <- T1 + T2 + T3

# True model parameters and change size
delta <- 0.3
theta1 <- c(-0.5, 0.6, 1) 
theta2 <- c(-0.5 + delta, 0.6 + delta, 1) 
ini <- (theta1 + theta2)/2

# Simulate data
y <- rbind(sim.y(theta=c(theta1[1:2], 0.5, theta1[3]), S.dist, T=T1, T.burn=100),
			     sim.y(theta=c(theta2[1:2], 0.5, theta2[3]), S.dist, T=T2, T.burn=100),
			     sim.y(theta=c(theta1[1:2], 0.5, theta1[3]), S.dist, T=T3, T.burn=100))


### Parameters for implementing CLDML
# Set spatial and temporal lag in computing CL
s.lag <- t.lag <- 1

# Calculate average number of times an obs used in CL for Ck
tmp <- D.cal(y, S.dist, s.lag, t.lag)
Ck <- mean(tmp$use.obs)
remedy <- tmp$remedy # number of marginal likelihood needed for correcting edge effect

# Calculate K for pruning step
p.min <- p.max <- 3
xi.min <- xi.max <- 0
K <- Ck*(log(S*T)*(p.min/2 - p.max) + (2 + p.max)*log(2) + xi.min - 2*xi.max - log(T))


### Run PELT for minimizing CLMDL
### One can run PELT with res=2 to further speed up computation
# (res=2 means cp can only happen on t such that t%%2==0, i.e. t is even.)
t0 <- proc.time()
pelt_result <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K, s.lag=s.lag,
                    t.lag=t.lag, min.length=0.1*T, ini=ini, res=2)
t_compute <- proc.time() - t0


# > t_compute
# user  system elapsed 
# 3039.74 15.94 4415.69 
# > pelt_result
# [1]   0  50 100 150

save.image("example_mcp.RData")
