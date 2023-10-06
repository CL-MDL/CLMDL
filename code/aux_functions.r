# Matern correlation
Matern <- function(d, rho, nu=0.5) {
  d1 <- sqrt(2*nu)*d/rho
  result <- 2^(1-nu)/gamma(nu)*d1^nu*besselK(d1, nu)
  result[d==0] <- 1
  result 
}


# Simulate ONE segment
sim.y <- function(theta, S.dist, T, T.burn=100) {
  Sigma <- theta[4]*Matern(S.dist, rho=theta[2], nu=theta[3])
  y <- matrix(0, nrow=T+T.burn, ncol(S.dist))
  e <- rmvnorm(n=T+T.burn, sigma=Sigma)
  y[1,] <- e[1,]
  for (i in 2:(T+T.burn)) {
    y[i,] <- theta[1]*y[i-1,] + e[i,]
  }
  y[-c(1:T.burn),]
}


# Calculate -logPL and D.pair with given distance matrix, lags and remedy matrix
pl <- function(theta, y, S.dist, s.lag=1, t.lag=1, remedy) {
  # All pairs combinations within s.lag distance
  # No truncation for time.lag > 0
  s.dist.pair <- function(S.dist, s.lag) {
    tmp <- cbind(c(row(S.dist)*(S.dist<=s.lag)), c(col(S.dist)*(S.dist<=s.lag)), c(S.dist))
    tmp[tmp[,1]!=0,]
  }
  # Truncate for time.lag = 0 (remove repetition)
  s.dist.pair.0 <- function(S.dist,s.lag) {
    S.dist2 <- row(S.dist) - col(S.dist) < 0
    tmp <- cbind(c(row(S.dist)*(S.dist<=s.lag)), c(col(S.dist)*(S.dist<=s.lag)), c(S.dist))[c(S.dist2),]
    tmp[tmp[,1]!=0,]
  }
  
  S <- ncol(y) ; T <- nrow(y)
  phi <- theta[1] ; rho <- theta[2] ; nu <- theta[3] ; sigma2 <- theta[4]
  logPL <- D.pair <- 0

  Sigma <- theta[4]*Matern(S.dist, rho=theta[2], nu=theta[3])
  y.var <- sigma2/(1-phi^2)
  
  s.comb <- s.dist.pair(S.dist, s.lag)
  s.comb.0 <- s.dist.pair.0(S.dist, s.lag)
  s.len <- unique(s.comb[,3])[-1] # s.comb includes 0, so remove 0 for non-zero spatial lag
  
  # Time lag 1 to max t lag
  for (t in 1:t.lag) {
  	# Space-lag 0
  	S.tmp <- matrix(y.var*c(1, phi^t, phi^t, 1), ncol=2)
  	for (j in 1:S){
    	x.tmp <- cbind(y[1:(T-t),j], y[(t+1):T,j])
    	logPL <- logPL + sum(dmvnorm(x=x.tmp, sigma=S.tmp, log=TRUE) + log(2*pi))
  	}
  	# Space-lag > 0
  	for (s in s.len) {
		  coord.tmp <- s.comb[s.comb[,3]==s,]
		  cov.tmp <- phi^t/(1-phi^2)*Sigma[coord.tmp[1,1], coord.tmp[1,2]]
		  S.tmp <- matrix(c(y.var, cov.tmp, cov.tmp, y.var), ncol=2)
		
  		for (i in 1:nrow(coord.tmp)) {
    		x.tmp <- cbind(y[1:(T-t), coord.tmp[i,1]], y[(t+1):T, coord.tmp[i,2]])
    		logPL <- logPL + sum(dmvnorm(x=x.tmp, sigma=S.tmp, log=TRUE) + log(2*pi))
  		}
  	}
  }
  # Time-lag 0, Space-lag > 0
  for (s in s.len) {
	  coord.tmp <- s.comb.0[s.comb.0[,3]==s,]
	  cov.tmp <- Sigma[coord.tmp[1,1], coord.tmp[1,2]]/(1-phi^2)
	  S.tmp <- matrix(c(y.var, cov.tmp, cov.tmp, y.var), ncol=2)
	
  	for (i in 1:nrow(coord.tmp)) {
    	x.tmp <- cbind(y[,coord.tmp[i,1]], y[,coord.tmp[i,2]])
    	logPL <- logPL + sum(dmvnorm(x=x.tmp, sigma=S.tmp, log=TRUE) + log(2*pi))
  	}
  }
  # Remedy for edge effect via marginal likelihood
  sd.remedy <- sqrt(y.var)
  for (i in 1:t.lag) {
    logPL <- logPL + sum(dnorm(rep(y[i,], c(remedy[i,])), sd=sd.remedy, log=TRUE) + log(2*pi)/2) + # remedy at start of segment
    				 sum(dnorm(rep(y[(T-i+1),], c(remedy[i,])), sd=sd.remedy, log=TRUE) + log(2*pi)/2)	# remedy at end of segment
  }
  list(nlogPL = -logPL)
}


pl.nu <- function(theta, y, S.dist, s.lag=1, t.lag=1, remedy) {
  theta <- c(theta[1:2], 0.5, theta[3])
  pl(theta, y, S.dist, s.lag, t.lag, remedy)
}


nlogPL.nu <- function(theta, y, S.dist, s.lag=1, t.lag=1, remedy) {
  # cat(theta,"\n")
  theta <- c(theta[1:2], 0.5, theta[3])
  pl(theta, y, S.dist, s.lag, t.lag, remedy)$nlogPL
}


# Input y and distance matrix find the D.pair (total number of pairs), use.obs (no. of time each obs used) and
# a remedy matrix for the edge effect
D.cal <- function(y, S.dist, s.lag=1, t.lag=1) {
  # All pairs combinations within s.lag distance
  # No truncation for time.lag > 0
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
  S <- ncol(y); T <- nrow(y)
  D.pair <- 0
  
  obs <- 0*y
  s.comb <- s.dist.pair(S.dist, s.lag)
  s.comb.0 <- s.dist.pair.0(S.dist, s.lag)
  s.len <- unique(s.comb[,3])[-1]
  
  # Time lag 1 to max t lag
  for (t in 1:t.lag) {
  	# Space-lag 0
  	for (j in 1:S){
    	obs[1:(T-t), j] <- obs[1:(T-t), j] + 1
    	obs[(t+1):T, j] <- obs[(t+1):T, j] + 1
    	x.tmp <- cbind(y[1:(T-t), j], y[(t+1):T, j])
    	D.pair <- D.pair + nrow(x.tmp)
  	}
  	# Space-lag > 0
  	for (s in s.len) {
  	  coord.tmp <- s.comb[s.comb[,3]==s,]
  		for (i in 1:nrow(coord.tmp)) {
    		obs[1:(T-t), coord.tmp[i,1]] <- obs[1:(T-t), coord.tmp[i,1]] + 1
    		obs[(t+1):T, coord.tmp[i,2]] <- obs[(t+1):T, coord.tmp[i,2]] + 1 
    		x.tmp <- cbind(y[1:(T-t), coord.tmp[i,1]], y[(t+1):T, coord.tmp[i,2]])
    		D.pair <- D.pair + nrow(x.tmp)
  		}
  	}
  }
  # Time-lag 0, Space-lag > 0
  for (s in s.len) {
  	coord.tmp <- s.comb.0[s.comb.0[,3]==s,]
  	for (i in 1:nrow(coord.tmp)) {
    	obs[,coord.tmp[i,1]] <- obs[,coord.tmp[i,1]] + 1
    	obs[,coord.tmp[i,2]] <- obs[,coord.tmp[i,2]] + 1
    	x.tmp <- cbind(y[,coord.tmp[i,1]],y[,coord.tmp[i,2]])
    	D.pair <- D.pair + nrow(x.tmp)
  	}
  }
  # By symmetric first row of remedy stands for the number of time of marginal for 1 and T, 
  # (if applicable) second row for 2 and T-1 and so on.
  remedy <- matrix(NA, nrow=t.lag, ncol=S)
  for (i in 1:t.lag) {
    remedy[i,] <- obs[t.lag+1,]-obs[i,]
  }
  D.pair <- D.pair + sum(remedy)
  list(D.pair=D.pair, use.obs=obs, remedy=remedy)
}


# Cost function for a segment
cost <- function(y.full, start, end, T, S.dist, s.lag=1, t.lag=1, Ck, min.length, remedy, ini) {
  S <- ncol(y.full)
  if (end - start + 1 < min.length) {
    return(Inf) 
  } else {
    y.seg <- y.full[start:end,]
    vec.par <- optim(ini, nlogPL.nu, y=y.seg, S.dist=S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy,
    		method="L-BFGS-B", lower=c(-0.7,0.1,1e-3), upper=c(0.7,3,5))$par
    
    obj <- pl.nu(theta=vec.par, y=y.seg, S.dist, s.lag=s.lag, t.lag=t.lag, remedy=remedy)
    lik <- obj$nlogPL
    pen <- Ck*(3/2*(log(end-start+1) + log(S)) + log(end-start+1))
  }
  lik + pen
}


# PELT function
pelt <- function(y.full, S.dist, Ck, remedy, K, s.lag=1, t.lag=1, min.length, ini, res=1) {
  T <- nrow(y.full)
  S <- ncol(y.full)

  # For initializing the objective function F, the cp set cp[0] and the candidate set R[0]
  F <- c(0, rep(0,T))
  cp <- 0
  R <- list(0)
  
  # For minimum segment length initialization
  for (i in 1:(2*min.length-1)) {
  	F[i+1] <- cost(y.full=y.full, start=1, end=i, T=T, S.dist=S.dist,
  	               s.lag=s.lag, t.lag=t.lag, Ck=Ck, min.length=min.length, remedy=remedy, ini=ini)
  	cpi <- list(c(0,i))
  	cp <- c(cp, cpi)
  	R.next <- list(c(0, i-min.length+1))
  	R <- c(R, R.next)
  	cat(i,"\n")
  }
  # Main function
  for(i in (2*min.length):T) {
  	possible.F <- rep(0, length(R[[i]]))
  	cnt <- 1
  	for(j in R[[i]]) {
    	possible.F[cnt] <- F[j+1] + cost(y.full=y.full, start=j+1, end=i, T=T, S.dist=S.dist,
    	                                 s.lag=s.lag, t.lag=t.lag, Ck=Ck, min.length=min.length,
    	                                 remedy=remedy, ini=ini)
    	cnt <- cnt+1
  	}
  	F[i+1] <- min(possible.F)
  	tau1 <- R[[i]][possible.F==F[i+1]][1]
  	cpi <- list(c(cp[[tau1+1]], i))
  	cp <- c(cp, cpi)
  	
  	cat(i, "est cp:", cp[[i+1]], " candidate:", R[[i]], "\n")	# The current "best" solution by considering the first i data-points
  	
  	# Pruning step in PELT
  	R.next <- c(R[[i]][possible.F+K<=F[i+1]], i-min.length+1)
  	# Further speed up computation by only considering every res time candidate
  	# By default res=1 and thus skip this speed up
  	if(res!=1){
  	  R.next <- R.next[R.next%%res==0]
  	}
  	R <- c(R, list(R.next))
  }
  cp[[T+1]] 						# Final output is the change-point position by considering the whole dataset.
}


# Given change point location, calculate CLMDL (cost function)
total_cost <- function(cp_vec, y.full, T, S.dist, s.lag=1, t.lag=1, Ck, remedy, min.length, ini) {
	ans <- 0
	len <- length(cp_vec)
	for (i in 1:(len-1)) {
		ans <- ans + cost(y.full, start=cp_vec[i]+1, end=cp_vec[i+1], 
				              T, S.dist, s.lag, t.lag, Ck, min.length, remedy, ini)
	}
	ans
}


# CI function
Q.grid <- function(lambda.hat, T, S.dist, y.full, theta1, theta2, s.lag, t.lag, remedy, q.bound) {
  q.range <- -q.bound:q.bound
  q.range <- q.range[q.range+lambda.hat >=1 & q.range+lambda.hat <= T]
  cl <- 0*q.range
  q_min <- max(1, min(q.range) + lambda.hat - 5)
  q_max <- min(T, max(q.range) + lambda.hat + 5)
  
  for (i in 1:length(cl)) {
    cp_tmp <- lambda.hat + q.range[i]
    cl[i] <- pl.nu(theta=theta1, y=y.full[q_min:cp_tmp,], S.dist, s.lag, t.lag, remedy)$nlogPL +
             pl.nu(theta=theta2, y=y.full[(cp_tmp+1):q_max,], S.dist, s.lag, t.lag, remedy)$nlogPL
  }
  q.range[cl==min(cl)]
}

