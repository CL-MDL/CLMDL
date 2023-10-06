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
T1 <- T2 <- 50
T <- T1 + T2

# True model parameters and change size
delta <- 0.3
theta1 <- c(-0.5, 0.6, 1)
theta2 <- c(-0.5 + delta, 0.6 + delta, 1) 
ini <- (theta1 + theta2)/2

# Simulate data
y <- rbind(sim.y(theta=c(theta1[1:2], 0.5, theta1[3]), S.dist=S.dist, T=T1, T.burn=100),
			     sim.y(theta=c(theta2[1:2], 0.5, theta2[3]), S.dist=S.dist, T=T2, T.burn=100))


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
t0 <- proc.time()
pelt_result <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K, s.lag=s.lag,
                    t.lag=t.lag, min.length=0.1*T, ini=ini, res=1)
t_compute <- proc.time() - t0

# > pelt_result
# [1]   0  50 100
# > t_compute
# user  system elapsed 
# 4071.11 4.45 4480.97 


### Subsequent analysis after CP detection
# first segment
y.seg1 <- y[(pelt_result[1]+1):pelt_result[2],]
par1 <- optim(ini, nlogPL.nu, y=y.seg1, S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy,
    		      method="L-BFGS-B", lower=c(-0.7,0.1,1e-3), upper=c(0.7,3,5))$par	
cost1 <- cost(y, start=pelt_result[1]+1, end=pelt_result[2], T=T, 
			        S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, Ck=Ck, min.length=0.1*T, remedy=remedy, ini=ini)
			
# second segment
y.seg2 <- y[(pelt_result[2]+1):pelt_result[3],]
par2 <- optim(ini, nlogPL.nu, y=y.seg2, S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy,
    		      method="L-BFGS-B", lower=c(-0.7,0.1,1e-3), upper=c(0.7,3,5))$par	
cost2 <- cost(y, start=pelt_result[2]+1, end=pelt_result[3], T=T, 
			        S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, Ck=Ck, min.length=0.1*T, remedy=remedy, ini=ini)

# check total cost
cost1 + cost2
total_cost(pelt_result, y.full=y, T=T, S.dist=S.dist, s.lag=s.lag, t.lag=t.lag,
           Ck=Ck, min.length=10, remedy=remedy, ini=ini)


### CI construction
B <- 100
theta1_est <- par1
theta2_est <- par2
lambda_est <- pelt_result[-c(1,length(pelt_result))] # estimated cp location
q_emp <- c()
t0 <- proc.time()
for(b in 1:B){
  set.seed(b)
  print(b)
  yB <- rbind(sim.y(theta=c(theta1_est[1:2], 0.5, theta1_est[3]), S.dist=S.dist, T=lambda_est, T.burn=100),
              sim.y(theta=c(theta2_est[1:2], 0.5, theta2_est[3]), S.dist=S.dist, T=T-lambda_est, T.burn=100))
  
  q_emp[b] <- Q.grid(lambda.hat=lambda_est, T=T, S.dist=S.dist, y.full=yB, theta1=theta1_est, theta2=theta2_est,
                     s.lag=s.lag, t.lag=t.lag, remedy=remedy, q.bound=round(0.2*T))
}
t_compute3 <- proc.time() - t0

ci <- quantile(q_emp, probs=c(0.025,0.975))

# > ci
# 2.5%  97.5% 
# -5.050  4.525 
# 
# > t_compute3
# user  system elapsed 
# 118.25    0.14  182.88 


### One can run PELT with res=2 to further speed up computation
# (res=2 means cp can only happen on t such that t%%2==0, i.e. t is even.)
t0 <- proc.time()
pelt_result2 <- pelt(y.full=y, S.dist=S.dist, Ck=Ck, remedy=remedy, K=K, s.lag=s.lag,
                     t.lag=t.lag, min.length=0.1*T, ini=ini, res=2)
t_compute2 <- proc.time() - t0

# > pelt_result2
# [1]   0  50 100
# > t_compute2
# user  system elapsed 
# 1672.68  11.52 2437.55


save.image('example_main.RData')
