rm(list=ls())
library(mvtnorm)
library(geosphere)
library(dplyr)
library(purrr)
library(pbv)
library(foreach)
library(doParallel)
source("./Util_RealData.R")
############# Code for real data analysis in the manuscript


### Number of cores used for parallel computing in pelt_parallel()
num_cores <- 100


### Read in data
load("./station_month_76.RData")
TT <- nrow(data.mat)
S <- ncol(data.mat)
S.dist <- dist.geosphere(info.select[,3:2])  # longitude, latitude
y_star <- as.matrix(data.mat)


### Preprocessing of the raw data
# log transformation
y.full <- log(y_star+1)

# stationarizing transformation for mu and sigma (see Lund et al 1995)
month_row <- rep(1:12, times=60)
mu_matrix <- sigma_matrix <- matrix(NA, nrow=12, ncol=76)

for (j in 1:76) {
	for (i in 1:12) {
		tmp <- y.full[month_row==i, j]
		mu_matrix[i,j] <- mean(tmp) # station-month mean
		sigma_matrix[i,j] <- sd(tmp) # station-month sd
	}
}

mu_matrix <- do.call(rbind, replicate(60, mu_matrix, simplify=FALSE))
sigma_matrix <- do.call(rbind, replicate(60, sigma_matrix, simplify=FALSE))
y.full <- (y.full-mu_matrix)/sigma_matrix


### Set up covariates in the linear regression for the mean function
lat <- matrix(rep(info.select$lat, each=TT), nrow=TT)
long <- matrix(rep(info.select$long, each=TT), nrow=TT)
elevation	<- matrix(rep(info.select$elev, each=TT), nrow=TT)
X.full <- list(lat, long, elevation)
ini <- c(3, 0.5, 1.5, 1) # initial point for CLMLE parameter estimation


### Parameters for implementing CLDML
# Set spatial and temporal lag in computing CL
s.lag <- 500
t.lag <- 3
min.length <- 150 # minimum length between two consecutive change-point

# Collect all pairs of observations within s.lag distance (will be used in PL)
s.comb <- data.frame(s.dist.pair(S.dist, s.lag))
names(s.comb) <- c("ind1", "ind2", "distance")
s.comb.0 <- s.dist.pair.0(S.dist, s.lag)
s.len <- unique(s.comb[,3])
s.len <- s.len[s.len>0] # all unique spatial distance

res <- D.cal(y=y.full, S.dist, s.lag, t.lag)
Ck <- mean(res$use.obs)
remedy <- res$remedy # number of marginal likelihood needed for correcting edge effect
p.min <- p.max <- length(ini)+1+length(X.full)
K <- Ck*(log(S*TT)*(p.min/2-p.max)+(2+p.max)*log(2)-log(TT)) # K for pruning step


######### Main function ##########
### 1. Estimate change-points based on PELT via minimizing CLMDL
### The PELT is written here with parallel computing
### to speed up the dynamic programming
t0 <- proc.time()
est_cp_result <- pelt_parallel(y.full=y.full, X.full=X.full, S.dist=S.dist,
                               Ck=Ck, remedy=remedy, K=K, s.lag=s.lag, t.lag=t.lag,
                               min.length=min.length, num_cores=num_cores)
t_compute_est <- proc.time() - t0
pelt_result <- est_cp_result$cp[[length(est_cp_result$cp)]] # estimated change-points
# > pelt_result
# [1]   0 150 327 720

# > t_compute_est[3]
# elapsed
# 15487.74 


### 2. Parameter estimation based on estimated change-points
num_seg <- length(pelt_result)-1
est_para_result <- list()
for(seg_index in 1:num_seg){
  start <- pelt_result[seg_index]+1
  end <- pelt_result[seg_index+1]
  
  y.seg <- y.full[start:end, ]
  X.seg <- do.call(cbind, lapply(X.full, function(x, start, end) {c(x[start:end,])}, start=start, end=end))
  est_seg <- para_est(y=y.seg, X=X.seg, ini=ini, S.dist, s.lag, t.lag, remedy)
  
  est_para_result[[seg_index]] <- est_seg
}


### 3. CI Construction for the estimated change-points
ci_result <- q_result <- c()
simu_seg_length <- min(diff(pelt_result)) # length for simulated segment
num_cp <- num_seg-1
t0 <- proc.time()
for(cp_index in 1:num_cp){
  prev_cp <- pelt_result[cp_index]
  cur_cp <- pelt_result[cp_index+1]
  next_cp <- pelt_result[cp_index+2]
  
  # Parameter estimation on the pre-change segment
  start1 <- max(cur_cp-simu_seg_length+1, prev_cp+1)
  end1 <- cur_cp
  y.seg1 <- y.full[start1:end1, ]
  X.seg1 <- do.call(cbind, lapply(X.full, function(x, start, end) {c(x[start:end,])}, start=start1, end=end1))
  tmp_est_seg1 <- para_est(y=y.seg1, X=X.seg1, ini=ini, S.dist, s.lag, t.lag, remedy)
  
  # Parameter estimation on the post-change segment
  start2 <- cur_cp+1
  end2 <- min(cur_cp+simu_seg_length, next_cp)
  y.seg2 <- y.full[start2:end2, ]
  X.seg2 <- do.call(cbind, lapply(X.full, function(x, start, end) {c(x[start:end,])}, start=start2, end=end2))
  tmp_est_seg2 <- para_est(y=y.seg2, X=X.seg2, ini=ini, S.dist, s.lag, t.lag, remedy)
  
  seg1_tmp <- nrow(y.seg1)
  seg2_tmp <- nrow(y.seg2)
  TT_tmp <- seg1_tmp + seg2_tmp
  
  # compute Cholesky decomposition for the non-separable Cressie-Huang covariance function
  beta1 <- tmp_est_seg1$beta
  theta1 <- tmp_est_seg1$theta
  beta2 <- tmp_est_seg2$beta
  theta2 <- tmp_est_seg2$theta
  chol_cov1 <-  chol_Cressie(theta1, S.dist, T_len=seg1_tmp) 
  chol_cov2 <-  chol_Cressie(theta2, S.dist, T_len=seg2_tmp)
  
  # Main function for CI construction
  B <- 100
  set.seed(12345)

  q_emp <- c()
  for (b in 1:B) {
    y.sim1 <- sim.y(X.seg1, beta1, theta1, S.dist, chol_cov=chol_cov1)
    y.sim2 <- sim.y(X.seg2, beta2, theta2, S.dist, chol_cov=chol_cov2)
    y.sim <- rbind(y.sim1, y.sim2)
    
    X <- do.call(cbind, lapply(X.full, function(x, start, end) {c(x[start:end,])}, start=start1, end=end2))
    q_emp[b] <- Q.grid(lambda.hat=seg1_tmp, TT=TT_tmp, S.dist=S.dist, X=X, y=y.sim, 
                       beta1=beta1, beta2=beta2, theta1=theta1, theta2=theta2, 
                       s.lag=s.lag, t.lag=t.lag, remedy=remedy, q.prop=0.3)
  }

  q_result <- cbind(q_result, q_emp)
  ci_lower_tmp <- -ceiling(-quantile(q_emp, probs=0.005))-1
  ci_upper_tmp <- ceiling(quantile(q_emp, probs=0.995))+1
  ci_result <- rbind(ci_result, pelt_result[cp_index+1]+c(ci_lower_tmp, ci_upper_tmp))
  
  rm(chol_cov1)
  rm(chol_cov2)
}
t_compute_ci <- proc.time() - t0

# > ci_result
# 0.5% 99.5%
# [1,]  136   170
# [2,]  303   348

# > t_compute_ci[3]
# elapsed 
# 6991.763 
