# Distance matrix by geosphere::distGeo
dist.geosphere <- function(longlat, round_para=NULL) {
  len <- nrow(longlat)
  ans <- matrix(0,ncol=len, nrow=len)
  for (i in 1:(len-1)) {
  	for (j in (i+1):len) { 
  		ans[i,j] <- distGeo(longlat[i,], longlat[j,])
  	}
  }
  ans <- (ans + t(ans))/10^3
  if (!is.null(round_para)) {
    ans <- round(ans/round_para)*round_para
  }
  ans
}


# Auxiliary function to compute pairs within s.lag distance
# No truncate for time.lag > 0
s.dist.pair <- function(S.dist, s.lag) {
  tmp <- cbind(c(row(S.dist)*(S.dist<=s.lag)), c(col(S.dist)*(S.dist<=s.lag)), c(S.dist))
  tmp[tmp[,1]!=0,]
}
# Truncate for time.lag = 0
s.dist.pair.0 <- function(S.dist, s.lag) {
  S.dist2 <- row(S.dist) - col(S.dist) < 0
  tmp <- cbind(c(row(S.dist)*(S.dist<=s.lag)), c(col(S.dist)*(S.dist<=s.lag)), c(S.dist))[c(S.dist2),]
  tmp[tmp[,1]!=0,]
}


# Input y and distance matrix find the D.pair (total number of pairs), use.obs (no. of time each obs used) and
# a remedy matrix for the edge effect
D.cal <- function(y, S.dist, s.lag, t.lag) {
  S <- ncol(y); TT <- nrow(y)
  D.pair <- 0
  
  obs <- 0*y
  s.comb <- s.dist.pair(S.dist, s.lag)
  s.comb.0 <- s.dist.pair.0(S.dist, s.lag)
  s.len <- unique(s.comb[,3])[-1]
  
  # Time lag 1 to max t lag
  for (t in 1:t.lag) {
    # Space-lag 0
    for (j in 1:S){
      obs[1:(TT-t),j] <- obs[1:(TT-t),j] + 1
      obs[(t+1):TT,j] <- obs[(t+1):TT,j] + 1
      x.tmp <- cbind(y[1:(TT-t),j], y[(t+1):TT,j])
      D.pair <- D.pair + nrow(x.tmp)
    }
    # Space-lag > 0
    for (s in s.len) {
      coord.tmp <- s.comb[s.comb[,3]==s,]
      for (i in 1:nrow(coord.tmp)) {
        obs[1:(TT-t),coord.tmp[i,1]] <- obs[1:(TT-t),coord.tmp[i,1]] + 1
        obs[(t+1):TT,coord.tmp[i,2]] <- obs[(t+1):TT,coord.tmp[i,2]] + 1 
        x.tmp <- cbind(y[1:(TT-t),coord.tmp[i,1]],y[(t+1):TT,coord.tmp[i,2]])
        D.pair <- D.pair + nrow(x.tmp)
      }
    }
  }
  # Time-lag 0, Space-lag > 0
  for (s in s.len) {
    coord.tmp <- s.comb.0[s.comb.0[,3] ==s,]
    if (NROW(coord.tmp) > 0 ){
      if (class(coord.tmp) == "numeric") {coord.tmp <- matrix(coord.tmp, ncol=3)}
      for (i in 1:NROW(coord.tmp)) {
        obs[,coord.tmp[i,1]] <- obs[,coord.tmp[i,1]] + 1
        obs[,coord.tmp[i,2]] <- obs[,coord.tmp[i,2]] + 1
        x.tmp <- cbind(y[,coord.tmp[i,1]],y[,coord.tmp[i,2]])
        D.pair <- D.pair + nrow(x.tmp)
      }
    }
  }
  # By symmetry, first row of remedy stands for the number of time of marginal for 1 and T, 
  # second row for 2 and T-1 and so on.
  remedy <- matrix(NA, nrow=t.lag, ncol=S)
  for (i in 1:t.lag) {
    remedy[i,] <- obs[t.lag+1,]-obs[i,]
  }
  D.pair <- D.pair + sum(remedy) 
  list(D.pair=D.pair, use.obs=obs, remedy=remedy)
}


# Cressie and Huang 1999 covariance function
# h: space distance
# u: time distance
# eta: (a, b, c, nu, sigma^2), all >(=) 0
Cressie <- function(h, u, eta) {
	a <- eta[1]; b <- eta[2]; c <- eta[3]
	nu <- eta[4]; sigma2 <- eta[5]

	au1 <- a^2*u^2 + 1
	auc <- a^2*u^2 + c
	bh <- b*sqrt(au1/auc)*h

  result <- sigma2*c/(au1^nu*auc)*exp(-bh/sqrt(2*nu))
	# result <-	sigma2*2*c/(au1^nu*auc*gamma(nu))*(bh/2)^nu*besselK(bh, nu) # same result when nu=0.5
	# result[h==0] <- (sigma2*c/(au1^nu*auc))[h==0]
	result
}


# log density of bivariate normal distribution
log_dbinorm <- function(z1, z2, sigma2, cov_z) {
  if(length(cov_z)!=NROW(z1)) { 
    rho <- rep(cov_z[1]/sigma2, NROW(z))
  } else {
    rho <- cov_z/sigma2
  }
  z1 <- z1/sqrt(sigma2)
  z2 <- z2/sqrt(sigma2)
  
  pbv::pbv_rcpp_dbvnorm(z1, z2, rho=rho, use_log=TRUE)-log(sigma2) # last term is Jacobian
}


# Calculate -logPL and D.pair with given distance matrix, lags and remedy matrix
pl <- function(theta, y, S.dist, s.lag=1, t.lag=1, remedy) {
  y <- as.matrix(y)
  S <- ncol(y); TT <- nrow(y)
  eta <- theta
  
  s.comb <- data.frame(s.dist.pair(S.dist, s.lag))
  names(s.comb) <- c("ind1", "ind2", "distance")
  s.comb.0 <- s.dist.pair.0(S.dist, s.lag)
  s.len <- unique(s.comb[,3])[-1] # s.comb includes 0, so remove 0 for non-zero spatial lag

  # First part of likelihood: Space-lag 0, time-lag 1,...,t.lag
  T_tmp <- diag(rep(1, TT))
  row_T <- c(row(T_tmp))
  col_T <- c(col(T_tmp))
  t_diff <- which(c(col_T-row_T) %in% 1:t.lag)
  tmp <- data.frame(t1=row_T[t_diff], t2=col_T[t_diff])
  cov_pair <- Cressie(h=0, u=tmp$t2-tmp$t1, eta)
  cov_pair <- rep(cov_pair, times=S)
  y1 <- c(y[tmp$t1,])
  y2 <- c(y[tmp$t2,])
  logPL <- sum(log_dbinorm(y1, y2, sigma2=eta[5], cov_pair))
  
  # Second part of likelihood: Space-lag > 0, time-lag 1,...,t.lag
  coord.within.s.lag <- s.comb %>% filter(distance<=s.lag, distance>0) %>% rename(h=distance)
  all_pair_combination <- purrr::map_dfr(seq_len(t.lag+1), ~coord.within.s.lag) %>%
	  mutate(u=rep(0:t.lag, each=NROW(coord.within.s.lag)))
	
  for(u_ind in 1:t.lag) {
  	tmp <- all_pair_combination %>% filter(u==u_ind) %>% mutate(cov_pair=Cressie(h, u, eta))
  	tmp2 <- data.frame(y1=c(y[1:(TT-u_ind), tmp$ind1]), y2=c(y[(u_ind+1):TT, tmp$ind2]),
  	                   cov_pair=rep(tmp$cov_pair, each=TT-u_ind)) %>%
  	          transmute(logPL=log_dbinorm(y1, y2, sigma2=eta[5], cov_pair)) %>% sum()
  	logPL <- logPL  + tmp2
  }
  
  # Third part of likelihood: Time-lag 0, Space-lag > 0
  tmp <- all_pair_combination %>% filter(u==0, ind1>ind2) %>% mutate(cov_pair=Cressie(h, u, eta))
  tmp3 <- data.frame(y1=c(y[,tmp$ind1]), y2=c(y[,tmp$ind2]), cov_pair=rep(tmp$cov_pair, each=TT)) %>% 
			      transmute(logPL=log_dbinorm(y1, y2, sigma2=eta[5], cov_pair)) %>% sum()
  logPL <- logPL + tmp3
  
  # Remedy for edge effect
  sd.remedy <- sqrt(eta[5])
  for (i in 1:t.lag) {
    logPL <- logPL + sum(dnorm(rep(y[i,], c(remedy[i,])), sd=sd.remedy, log=TRUE)) + # remedy at start of segment
      sum(dnorm(rep(y[(TT-i+1),], c(remedy[i,])), sd=sd.remedy, log=TRUE)) # remedy at end of segment
  }
  list(nlogPL=-logPL)
}


# Wrapper for pl()
nlogPL.fixed.nu <- function(theta, nu=0.5, y, S.dist, s.lag, t.lag, remedy=remedy) {
  if(any(theta < 1e-5)) ans <- Inf
  theta <- c(theta[1:3], nu, theta[-c(1:3)])
  ans <- pl(theta, y, S.dist, s.lag, t.lag, remedy=remedy)$nlogPL
  ans
}


# Cost function for a segment
cost <- function(y.full, X.full, start, end, TT, S.dist, s.lag=1, t.lag=1, Ck, min.length=150, remedy, ini) {
  S <- ncol(y.full) 
  if (end-start+1<min.length) {
    return(Inf) 
  } else {
    y.seg <- y.full[start:end,]
  	f.tmp <- function(x, start, end) {
  	  c(x[start:end,])
  	}
	  X.seg <- do.call(cbind, lapply(X.full, f.tmp, start=start, end=end))
	  names(X.seg) <- names(X.full)
		
	  # Stage 1: linear regression for the mean function
	  m.lm <- lm(c(y.seg)~X.seg)
    y.resid <- matrix(m.lm$resid, nrow=end-start+1, ncol=S)
	
    # Stage 2: Cressie-Huang covariance function estimation via CL
    obj <- optim(ini, nlogPL.fixed.nu, y=y.resid, nu=0.5, S.dist=S.dist, s.lag=s.lag,
                 t.lag=t.lag, remedy=remedy, method="L-BFGS-B", lower=rep(1e-4,4))
	  lik <- obj$value
	
	  # Penalty
	  num_para <- length(ini)+length(m.lm$coef)
    pen <- Ck*(num_para/2*(log(end-start+1)+log(S))+log(TT))
  }
  lik + pen
}


# PELT function
pelt_parallel <- function(y.full, X.full, S.dist, Ck, remedy, K, s.lag, t.lag, min.length,
                          ini=c(1, 1, 1, 1), num_cores=4) {
  registerDoParallel(num_cores)
  TT <- nrow(y.full)
  S <- ncol(y.full)
  
  # Initialization of PELT
  f <- c(0, rep(Inf, TT))
  cp <- 0
  R <- list(0)
  
  # For minimum segment length initialization
  init_cost <- foreach(i=1:(2*min.length-1), .combine=c, .packages=c("dplyr", "purrr", "pbv", "mvtnorm")) %dopar% {
    cost(y.full, X.full, start=1, end=i, TT, S.dist, s.lag, t.lag, Ck, min.length, remedy, ini)
  }
  for (i in 1:(2*min.length-1)) {
    f[i+1] <- init_cost[i]
    cpi <- list(c(0, i))
    cp <- c(cp, cpi)
    R.next <- list(c(0, i-min.length+1))
    R <- c(R, R.next)
  }
  
  # Main PELT (Pruned dynamic programming) for minimizing CLMDL
  for(i in (2*min.length):TT) {
    possible.F <- foreach(j=R[[i]], .combine=c, .packages=c("dplyr", "purrr", "pbv", "mvtnorm")) %dopar% {
      f[j+1] + cost(y.full, X.full, start=j+1, end=i, TT, S.dist, s.lag, t.lag, Ck, min.length, remedy, ini)
    }
    f[i+1] <- min(possible.F)
    tau1 <- R[[i]][possible.F==f[i+1]][1]
    cpi <- list(c(cp[[tau1+1]], i))
    cp <- c(cp, cpi)
    
    R.next <- list(c(R[[i]][possible.F+K<=f[i+1]], i-min.length+1)) # Pruned here
    
    R <- c(R, R.next)
  }
  result <- list(cp=cp, R=R)
}


# Parameter estimation on a given segment
para_est <- function(y, X, ini, S.dist, s.lag, t.lag, remedy){
  m.lm <- lm(c(y)~X)
  beta <- m.lm$coef # regression for the mean function
  y.resid <- matrix(m.lm$resid, nrow=nrow(y), ncol=ncol(y))
  theta <- optim(ini, nlogPL.fixed.nu, y=y.resid, 
                 nu=0.5, S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy,
                 method="L-BFGS-B", lower=rep(1e-4,4))$par # Cressie covariance function
  return(list(beta=beta, theta=theta))
}


# Simulate data via the Linear Regression mean function + Cressie and Huang (1995) covariance function
sim.y <- function(X.seg, beta, theta, S.dist, chol_cov=NULL) {
  S <- nrow(S.dist)
  ST <- nrow(X.seg)
  T_len <- ST/S
  # Mean function
  Xbeta <- matrix(c(cbind(1, X.seg)%*%matrix(beta, nc=1)), nc=S)
  
  if(is.null(chol_cov)) {
    chol_cov <- chol_Cressie(theta, S.dist, T_len)
  }
  # Plus noise with the Cressie-Huang covariance function
  Xbeta + matrix(c(chol_cov%*%rnorm(ST)), nc=S)
}


# Compute Cholesky decomposition for Cressie-Huang covariance function
chol_Cressie <- function(theta, S.dist, T_len) {
  S <- nrow(S.dist)
  eta <- c(theta[1:3], 0.5, theta[-c(1:3)]) # nu is fixed at 0.5
  h <- kronecker(S.dist, matrix(1, nr=T_len, nc=T_len))
  tmp <- matrix(1, nr=T_len, nc=T_len)
  tmp <- abs(row(tmp)-col(tmp))
  u <- kronecker(matrix(1, nr=S, nc=S), tmp)
  
  t(chol(Cressie(h, u, eta) + diag(1e-10, nrow(u))))
}


# Compute CI
Q.grid <- function(lambda.hat, TT, S.dist, X, y, beta1, beta2, theta1, theta2, s.lag, t.lag, remedy, q.prop=0.3) {
  S <- nrow(S.dist)
  q.range <- floor(-q.prop*TT):ceiling(q.prop*TT)
  q.range <- q.range[q.range+lambda.hat >=1 & q.range+lambda.hat <= TT]
  cl <- 0*q.range
  
  y.resid1 <- y - matrix(c(cbind(1,X)%*%matrix(beta1, nc=1)), nc=S)
  y.resid2 <- y - matrix(c(cbind(1,X)%*%matrix(beta2, nc=1)), nc=S)
  
  q_min <- max(1, min(q.range) + lambda.hat - 10)
  q_max <- min(TT, max(q.range) + lambda.hat + 10)
  
  for (i in 1:length(cl)) {
    cp_tmp <- lambda.hat + q.range[i]
    cl[i] <- nlogPL.fixed.nu(theta=theta1, nu=0.5, y=y.resid1[q_min:cp_tmp,], S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy)+
      nlogPL.fixed.nu(theta=theta2, nu=0.5, y=y.resid2[(cp_tmp+1):q_max,], S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy)
  }
  q.range[cl==min(cl)]
}

