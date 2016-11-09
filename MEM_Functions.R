library(BRugs)
library(gtools)
library(matrixStats)

#############################################
### General functions
#############################################

boa.hpd <- function(x, alpha){
###Calculate HPD Intervals given vector of values and alpha
    n <- length(x)
    m <- max(1, ceiling(alpha * n))
    y <- sort(x)
    a <- y[1:m]
    b <- y[(n - m + 1):n]
    i <- order(b - a)[1]
    structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}


#############################################
### Meta Analysis Approach Functions
#############################################
shm_calc_sim = function(x,S,N,scen,tau.up=50){
###calculate meta-analysis approach for estimation
###x: vector mean values for studies
###S: vector of standard deviations for studies
###N: sample size for studies
###scen: text string used to write out file for easy identification
###tau.up: upper limit of uniform prior on tau (used to calculate precision)

	mu.seq = seq(from=-15,to=10,length.out=500)
	asd = x[1]
	mu = mu.seq[asd]
	
	est.mu <- est.var <- est.mse <- est.bias <- est.esss <- est.esss.alt <- NULL #initialize vectors to store results
	est.xbar <- cov.ind <- hpd.width <- NULL

	for(k in 1:1000){
		set.seed(515+k)
		xbar = rnorm(1,mean=mu,sd=S[1]/sqrt(N[1])) #xbar.seq[k] #CENIC G
		xbar01 = x[2] #CENIC B
		xbar02 = x[3] #Quest Jr.
		xbar03 = x[4] #Quest Sr.

		s = S[1]^2 #variance
		s01 = S[2]^2
		s02 = S[3]^2
		s03 = S[4]^2

		n = N[1] #sample sizes
		n01 = N[2]
		n02 = N[3]
		n03 = N[4]

		v = s/n
		v01 = s01/n01
		v02 = s02/n02
		v03 = s03/n03
	
		data <- list(k=4,upper.unif=tau.up,y=c(xbar,xbar01,xbar02,xbar03),prec.y=c(1/v,1/v01,1/v02,1/v03))

		modelCheck("re_model_MetaAnalysis.txt")
		modelData(bugsData(data))
		modelCompile(numChains=3)
		modelInits('inits_MetaAnalysis.txt', chainNum=1)
		modelInits('inits_MetaAnalysis.txt', chainNum=2)
		modelInits('inits_MetaAnalysis.txt', chainNum=3)

		BURN <- 5000
		UPDATE <- 10000
		gC.fit <- BRugsFit(modelFile="re_model_MetaAnalysis.txt",data=data,numChains=3,para=c("theta"),nBurnin=BURN,nIter=UPDATE,DIC = FALSE, BRugsVerbose=FALSE)

		theta1 <- samplesSample('theta[1]')
		prec.theta1 <- 1/sd(theta1)^2

		mu.est <- mean(theta1)
		sd.comp = sd(theta1)
		bias.orig <- mu.est-mu
		esss.orig = n*(prec.theta1/(1/v)-1)

		hpd.bma_rep = boa.hpd(x=theta1,alpha=.05)

		include.mu = (hpd.bma_rep[1] < mu) & (mu < hpd.bma_rep[2])
		diff.hpd = hpd.bma_rep[2] - hpd.bma_rep[1]

		est.mu <- c(est.mu, mu.est)
		est.bias <- c(est.bias, bias.orig)
		est.esss <- c(est.esss, esss.orig)
		cov.ind = c(cov.ind, include.mu)
		hpd.width = c(hpd.width, diff.hpd)
	}

	var.part = var(est.mu)
	mean.part = mean(est.mu)
	bias.part = mean.part - mu
	mse.part = bias.part^2 + var.part
	esss.part = median(est.esss)
	esss.alt = s/var.part	
	abserr.part = abs(bias.part)
	cov.part = mean(cov.ind)
	hpd.part = mean(hpd.width)	

	mat = cbind(mu,mean.part,var.part,bias.part,mse.part,esss.part,esss.alt,abserr.part,cov.part,hpd.part)

	write.table(mat, file = paste0('ma.estimators_',scen,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)

#	return(mat)
}

shm_calc = function(x,S,N,tau.up=50){
###Function to calculate SHM results given data
###calculate meta-analysis approach for estimation
###x: vector mean values for studies
###S: vector of standard deviations for studies
###N: sample size for studies
###tau.up: upper limit of uniform prior on tau (used to calculate precision)

	set.seed(515)	

if(length(S)==4){
	xbar = x[1] #xbar.seq[k] #CENIC G
	xbar01 = x[2] #CENIC B
	xbar02 = x[3] #Quest Jr.
	xbar03 = x[4] #Quest Sr.

	s = S[1]^2 #variance
	s01 = S[2]^2
	s02 = S[3]^2
	s03 = S[4]^2

	n = N[1] #sample sizes
	n01 = N[2]
	n02 = N[3]
	n03 = N[4]

	v = s/n
	v01 = s01/n01
	v02 = s02/n02
	v03 = s03/n03
	
	data <- list(k=4,upper.unif=tau.up,y=c(xbar,xbar01,xbar02,xbar03),prec.y=c(1/v,1/v01,1/v02,1/v03))

	modelCheck("re_model_MetaAnalysis.txt")
	modelData(bugsData(data))
	modelCompile(numChains=3)
	modelInits('inits_MetaAnalysis.txt', chainNum=1)
	modelInits('inits_MetaAnalysis.txt', chainNum=2)
	modelInits('inits_MetaAnalysis.txt', chainNum=3)

	BURN <- 5000
	UPDATE <- 10000
	gC.fit <- BRugsFit(modelFile="re_model_MetaAnalysis.txt",data=data,numChains=3,para=c("theta"),nBurnin=BURN,nIter=UPDATE,DIC = FALSE, BRugsVerbose=FALSE)

	theta1 <- samplesSample('theta[1]')
	prec.theta1 <- 1/sd(theta1)^2

	mean.est <- mean(theta1)
	var.comp = var(theta1)
	esss.est = n*(prec.theta1/(1/v)-1)
}
if(length(S)==2){
	xbar = x[1] #xbar.seq[k] #CENIC G
	xbar01 = x[2] #CENIC B

	s = S[1]^2 #variance
	s01 = S[2]^2

	n = N[1] #sample sizes
	n01 = N[2]

	v = s/n
	v01 = s01/n01
	
	data <- list(k=2,upper.unif=tau.up,y=c(xbar,xbar01),prec.y=c(1/v,1/v01))

	modelCheck("re_model_MetaAnalysis.txt")
	modelData(bugsData(data))
	modelCompile(numChains=3)
	modelInits('inits_MetaAnalysis_2source.txt', chainNum=1)
	modelInits('inits_MetaAnalysis_2source.txt', chainNum=2)
	modelInits('inits_MetaAnalysis_2source.txt', chainNum=3)

	BURN <- 5000
	UPDATE <- 10000
	gC.fit <- BRugsFit(modelFile="re_model_MetaAnalysis.txt",data=data,numChains=3,para=c("theta"),nBurnin=BURN,nIter=UPDATE,DIC = FALSE, BRugsVerbose=FALSE)

	theta1 <- samplesSample('theta[1]')
	prec.theta1 <- 1/sd(theta1)^2

	mean.est <- mean(theta1)
	var.comp = var(theta1)
	esss.est = n*(prec.theta1/(1/v)-1)
}
	est.vec = c(mean.est, var.comp, esss.est)

	return(est.vec)
}


#############################################
### Commensurate Prior Approach Functions
#############################################

## No Borrowing variance for lambda: Hobbs, Sargent, Carlin, BA (12)
V.0 <- function(n,n1,sig2){ return( sig2/( n1*(1-(n1/n)) ) ) }  ## Equal to Var(yd.hat-y.hat)/( (1-n.d/n)^2 ) (12) & (15) ##

## Posterior Variance of lambda|tau ## (9)
V <- function(tau,n,n1,sig2,v0){ return( 1/( 1/(((n/n1)^2)*((sig2/n)+v0+(1/tau))) + ( ( n1*(1-(n1/n)) )/sig2 ) ) ) }

## Posterior mean for lambda|tau ## (10)
M <- function(tau,A,B,n,n1,sig2,v0){ return(  ( A/( (n/n1)*((sig2/n)+v0+(1/tau)) ) + ( B/(sig2/n1)) )/( 1/V(tau,n,n1,sig2,v0) ) ) }

## Shrinkage Weight ##
SW <- function(v,v0,tau){ return( ( v+v0 )/( v+v0+(1/tau) ) ) }

## Inverse Shrinkage Weight ##
invSW <- function(v,v0,w){ return( w/( (v+v0)*(1-w) ) ) }

## Hier. Precision from EHSS ##
#HierPrec <- function(EHSS,sig2,v0){ return(max(0.0000001,EHSS/( sig2-(v0*EHSS) ) )) }

HierPrec <- function(EHSS,sig2,v0){ 
#cat(EHSS, " ",sig2," ", v0, "\n")
 out <- EHSS/( sig2-(v0*EHSS) )
cat(out, "\n")
if(out<0){ EHSS. <- sig2/v0-0.0001
            out <- EHSS./( sig2-(v0*EHSS.) ) }
return(out) 
}

## Estimating tau for posterior ##
tau_est <- function(dat.C, dat.H, R, EmpBayes=0, Static=1, L.SLAB=0.01, SPIKE =10.0, U.SLAB=3, PI0=0.7, EHSS=1){
#####################################################################################
### Function for computing estimate for tau to use in posterior calculations ##################################
# Arguments:
# - dat.C = list of current data:
#		- n = total sample size
#		- nd = number assigned treatment, set to 0 for this paper
#		- ybar = overall sample mean
#		- yd.bar = sample mean of treatment
#		- s2 = model variance
# - dat.H = list of historical data:
#		- n = vector of sample sizes (one for each historical study)
#		- ybar = vector of sample means (one for each historical study)
#		- s2 = model variance for each historic source
# - R = number of future patients that remain to be randomized (used in AR.probability calculation)
# - EmpBayes = 0 - 1 indicator of Empirical Bayesian Inference
# - Static = 0 - 1 indicator of Static Borrowing Inference
# - U.SLAB = Hyperparameter for Spike-and-Slab Inference (default = 3)
# - PI0 = Hyperparameter for Spike-and-Slab Inference (default = 0.7)
# - W = shrinkage weight used for BOTH upper bound of MMLE tau in Empirical Bayes inference AND Shrinkage weight in Static Borrowing Inference (default = 0.5) ###
#

 ####################
 ## Data Constants ##
 ####################
 H <- length(dat.H$n)
 mu.hat <- (dat.C$n*dat.C$ybar-dat.C$nd*dat.C$yd.bar)/(dat.C$n-dat.C$nd)
 B <- dat.C$yd.bar-dat.C$ybar

 ##################################################################
 ## Fit to Concurrent Data alone corrected Hobbs et al 2012 (16) ##
 ##################################################################
  T.mu <- (dat.C$yd.bar-dat.C$ybar)/(1-dat.C$nd/dat.C$n)
  T.nu <- dat.C$n-2
  T.sig2 <- ( dat.C$s2 - (1/((1/dat.C$n)-(1/dat.C$nd)))*((dat.C$ybar-dat.C$yd.bar)^2) )/( dat.C$nd*T.nu*(1-dat.C$nd/dat.C$n) )
  VAR0.lambda <- T.sig2*(T.nu/(T.nu-2))

 ###########################
 ## MCMC for Joint Models ##
 ## 1) Spike & Slab       ##
 ## 2) Emp Bayes          ##
 ## 3) Static Borrowing   ##
 ###########################
 N.iter <- 15000; acc <- 0; P0 <- 0.5

 ## Fix sample variances at MMLEs (or REMLs) ##
 sig20 <- dat.H$s2 / dat.H$n ###ifelse(0.5*(dat.H$n-3)>1, dat.H$s2/((dat.H$n-3)-1), dat.H$s2/(dat.H$n-1) )
 sig2 <- dat.C$s2 / dat.C$n ###dat.C$s2/(dat.C$n-1)

 ## Set variance-dependent components ##
 omega0 <- sig20/dat.H$n
 v0 <- 1/sum( 1/omega0 )
 mu0.hat <- v0*sum( dat.H$ybar/omega0 )
 A <- dat.C$ybar - mu0.hat
 Delta.hat <- mu.hat - mu0.hat

 #####################
 ## Hyperparameters ##
 #####################
 ## Compute upper bound for tau from provided shrinkage weight or EHSS ##
 W <- SW( sig2/(dat.C$n-dat.C$nd), v0, HierPrec(EHSS,sig2,v0) )
 hyperParameters <- list(l.slab = HierPrec(L.SLAB,sig2,v0), 
                        u.slab = HierPrec(U.SLAB,sig2,v0), 
                        spike =  HierPrec(SPIKE,sig2,v0), 
                        pi0 = PI0, l.eb=0.005, u.eb=invSW( sig2/(dat.C$n-dat.C$nd), v0, W))

 ## Initial values (spike & slab) for mu, tau, & lambda ##
 tau <- hyperParameters$u.slab  
 lambda <- M(tau,A,B,dat.C$n,dat.C$nd,sig2,v0)
 mu <- ( (mu.hat/( sig2/(dat.C$n-dat.C$nd) )) + (mu0.hat/(v0+(1/tau))) )/( ((dat.C$n-dat.C$nd)/sig2) + (1/(v0+(1/tau))) )

 ## MMLE of tau, Hobbs et al 2012 (5) ##
 tau.EB <- 1/max( min( (Delta.hat^2)-(sig2/(dat.C$n-dat.C$nd))-v0, 1/hyperParameters$l.eb ), 1/hyperParameters$u.eb )

 return(tau.EB)
}


commensurateprior_calc_sim = function(x,S,N,scen){
###Function to generate simulation for commensurate prior
#x: vector mean values for studies
#S: vector of standard deviations for studies
#N: sample size for studies
#scen: text string used to write out file for easy identification

	mu.seq = seq(from=-15,to=10,length.out=500)
	asd = x[1]
	mu = mu.seq[asd]
	
	est.mu <- est.var <- est.mse <- est.bias <- est.esss <- est.esss.alt <- NULL #initialize vectors to store results
	est.xbar <- cov.ind <- hpd.width <- NULL

	for(k in 1:10000){
		set.seed(515+k)
		xbar = rnorm(1,mean=mu,sd=S[1]/sqrt(N[1])) #xbar.seq[k] #CENIC G
		xbar01 = x[2] #CENIC B
		xbar02 = x[3] #Quest Jr.
		xbar03 = x[4] #Quest Sr.

		s = S[1]^2 #variance
		s01 = S[2]^2
		s02 = S[3]^2
		s03 = S[4]^2

		n = N[1] #sample sizes
		n01 = N[2]
		n02 = N[3]
		n03 = N[4]

		v = s/n
		v01 = s01/n01
		v02 = s02/n02
		v03 = s03/n03

		dat.C <- list(n=n, nd=0, ybar=xbar, yd.bar=0, s2=s) #current study information
		dat.H<-list(n=c(n01,n02,n03),ybar=c(xbar01,xbar02,xbar03),s2=c(s01,s02,s03)) #historic sources information
		tau <- tau_est(dat.C=dat.C, dat.H=dat.H, R=0, EmpBayes=1, Static=0)
		tauinv <- 1/tau

		v0 <- (1/v01 + 1/v02 + 1/v03)^(-1)
		muhat.0 <- v0 * (xbar01/v01 + xbar02/v02 + xbar03/v03)

		mu.est <- xbar*((v0 + tauinv)/(v + v0 + tauinv)) + muhat.0*(v/(v + v0 + tauinv))
		var.comp <- (v0 + tauinv)*(v/(v + v0 + tauinv))
		bias.orig = (mu.est-mu)
		esss.orig = s/var.comp - n

		post.mu = rnorm(30000,mu.est,sd=sqrt(var.comp))
		hpd.bma_rep = boa.hpd(x=post.mu,alpha=.05)

		include.mu = (hpd.bma_rep[1] < mu) & (mu < hpd.bma_rep[2])
		diff.hpd = hpd.bma_rep[2] - hpd.bma_rep[1]

		est.mu <- c(est.mu, mu.est)
		est.bias <- c(est.bias, bias.orig)
		est.esss <- c(est.esss, esss.orig)
		cov.ind = c(cov.ind, include.mu)
		hpd.width = c(hpd.width, diff.hpd)
	}

	var.part = var(est.mu)
	mean.part = mean(est.mu)
	bias.part = mean.part - mu
	mse.part = bias.part^2 + var.part
	esss.part = median(est.esss) 
	esss.alt = s/var.part	
	abserr.part = abs(bias.part)
	cov.part = mean(cov.ind)
	hpd.part = mean(hpd.width)	

	mat = cbind(mu,mean.part,var.part,bias.part,mse.part,esss.part,esss.alt,abserr.part,cov.part,hpd.part)

	write.table(mat, file = paste0('commprior.estimators_',scen,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)
#	return(mat)
}

commensurateprior_calc = function(x,S,N){
###Function to calculate CP results given data
#x: vector mean values for studies
#S: vector of standard deviations for studies
#N: sample size for studies
	
if(length(S)==4){

	xbar = x[1]
	xbar01 = x[2] #CENIC B
	xbar02 = x[3] #Quest Jr.
	xbar03 = x[4] #Quest Sr.

	s = S[1]^2 #variance
	s01 = S[2]^2
	s02 = S[3]^2
	s03 = S[4]^2

	n = N[1] #sample sizes
	n01 = N[2]
	n02 = N[3]
	n03 = N[4]

	v = s/n
	v01 = s01/n01
	v02 = s02/n02
	v03 = s03/n03

	dat.C <- list(n=n, nd=0, ybar=xbar, yd.bar=0, s2=s) #current study information
	dat.H<-list(n=c(n01,n02,n03),ybar=c(xbar01,xbar02,xbar03),s2=c(s01,s02,s03)) #historic sources information
	tau <- tau_est(dat.C=dat.C, dat.H=dat.H, R=0, EmpBayes=1, Static=0)
	tauinv <- 1/tau

	v0 <- (1/v01 + 1/v02 + 1/v03)^(-1)
	muhat.0 <- v0 * (xbar01/v01 + xbar02/v02 + xbar03/v03)

	mean.est <- xbar*((v0 + tauinv)/(v + v0 + tauinv)) + muhat.0*(v/(v + v0 + tauinv))
	var.comp <- (v0 + tauinv)*(v/(v + v0 + tauinv))
	esss.est = s/var.comp - n
}
if(length(S)==2){

	xbar = x[1]
	xbar01 = x[2] 

	s = S[1]^2 #variance
	s01 = S[2]^2

	n = N[1] #sample sizes
	n01 = N[2]

	v = s/n
	v01 = s01/n01

	dat.C <- list(n=n, nd=0, ybar=xbar, yd.bar=0, s2=s) #current study information
	dat.H<-list(n=c(n01),ybar=c(xbar01),s2=c(s01)) #historic sources information
	tau <- tau_est(dat.C=dat.C, dat.H=dat.H, R=0, EmpBayes=1, Static=0)
	tauinv <- 1/tau

	v0 <- (1/v01)^(-1)
	muhat.0 <- v0 * (xbar01/v01)

	mean.est <- xbar*((v0 + tauinv)/(v + v0 + tauinv)) + muhat.0*(v/(v + v0 + tauinv))
	var.comp <- (v0 + tauinv)*(v/(v + v0 + tauinv))
	esss.est = s/var.comp - n
}
	est.vec = c(mean.est, var.comp, esss.est)

	return(est.vec)
}


#############################################
### MEM functions
#############################################

MEM_sim_calc = function(x,S,N,prior,scen,type){
###x: vector mean values for primary and then supplemental cohorts
#S: vector of standard deviations for primary and then supplemental cohorts
#N: sample size for primary and then supplemental cohorts
#prior: prior to use for calculation
#scen: text string used to write out file for easy identification

	mu.seq = seq(from=-15,to=10,length.out=500)
	asd = x[1]
	mu = mu.seq[asd]
	
	est.mu <- est.var <- est.mse <- est.bias <- est.esss <- est.esss.alt <- NULL #initialize vectors to store results
	est.xbar <- cov.ind <- hpd.width <- NULL

	for(k in 1:10000){
		set.seed(515+k)
		xbar = rnorm(1,mean=mu,sd=S[1]/sqrt(N[1])) #xbar.seq[k] #CENIC G
		xbar01 = x[2] #CENIC B
		xbar02 = x[3] #Quest Jr.
		xbar03 = x[4] #Quest Sr.

		s = S[1]^2 #variance
		s01 = S[2]^2
		s02 = S[3]^2
		s03 = S[4]^2

		n = N[1] #sample sizes
		n01 = N[2]
		n02 = N[3]
		n03 = N[4]

		v = s/n
		v01 = s01/n01
		v02 = s02/n02
		v03 = s03/n03
	
		###Posterior components
		##Means for each model
		M1 = xbar
		M2 = (v01*xbar + v*xbar01)/(v+v01)
		M3 = (v02*xbar + v*xbar02)/(v+v02)
		M4 = (v03*xbar + v*xbar03)/(v+v03)
		M5 = (v01*v02*xbar + v*v02*xbar01 + v*v01*xbar02)/(v01*v02 + v*v02 + v*v01)
		M6 = (v01*v03*xbar + v*v03*xbar01 + v*v01*xbar03)/(v01*v03 + v*v03 + v*v01)
		M7 = (v02*v03*xbar + v*v03*xbar02 + v*v02*xbar03)/(v02*v03 + v*v03 + v*v02)
		M8 = (v01*v02*v03*xbar + v*v02*v03*xbar01 + v*v01*v03*xbar02 + v*v01*v02*xbar03)/(v01*v02*v03 + v*v02*v03 + v*v01*v03 + v*v01*v02)

		##Variances for each model
		V1 = v
		V2 = (1/v + 1/v01)^-1
		V3 = (1/v + 1/v02)^-1
		V4 = (1/v + 1/v03)^-1
		V5 = (1/v + 1/v01 + 1/v02)^-1
		V6 = (1/v + 1/v01 + 1/v03)^-1
		V7 = (1/v + 1/v02 + 1/v03)^-1
		V8 = (1/v + 1/v01 + 1/v02 + 1/v03)^-1

		##Effective historical sample size for each model
		e1 = (s/V1) - n
		e2 = (s/V2) - n
		e3 = (s/V3) - n
		e4 = (s/V4) - n
		e5 = (s/V5) - n
		e6 = (s/V6) - n
		e7 = (s/V7) - n
		e8 = (s/V8) - n

		w <- calc.weights_MEM(xvec=c(xbar,x[2],x[3],x[4]),svec=c(S[1],S[2],S[3],S[4]),nvec=c(N[1],N[2],N[3],N[4]),prior=prior)
		w1<-w[1]; w2<-w[2]; w3<-w[3]; w4<-w[4]; w5<-w[5]; w6<-w[6]; w7<-w[7]; w8<-w[8]

		###Calculate MSE
		##Variance component
		#Weight denominator needed for delta method g function calculation (weight for each model):
		##Weight is w1 = 1/a1, calculations as a1 for ease of inclusion with quotient rule for derivatives
		a1 = (w1+w2+w3+w4+w5+w6+w7+w8)/w1 #note that 1/a1 equals the weight for model 1, and 1/a2 for model 2, etc.
		a2 = (w1+w2+w3+w4+w5+w6+w7+w8)/w2
		a3 = (w1+w2+w3+w4+w5+w6+w7+w8)/w3
		a4 = (w1+w2+w3+w4+w5+w6+w7+w8)/w4
		a5 = (w1+w2+w3+w4+w5+w6+w7+w8)/w5
		a6 = (w1+w2+w3+w4+w5+w6+w7+w8)/w6
		a7 = (w1+w2+w3+w4+w5+w6+w7+w8)/w7
		a8 = (w1+w2+w3+w4+w5+w6+w7+w8)/w8

		#Derivative of exponential part of each w_i weight component for delta method variance calculation
		##Note that the exp() part is left off because it is contained within the w_i part included in the next step
		##(i.e., (w1+w2+w3+...+w8)/w2 will contain all the necessary info besides the deriv of the exp part calculated for b here)
		b1 = 0
		b2 = -(xbar - xbar01)/(v + v01)
		b3 = -(xbar - xbar02)/(v + v02)
		b4 = -(xbar - xbar03)/(v + v03)
		b5 = -(xbar - xbar01)/(v + v01 + v*v01/v02) - (xbar - xbar02)/(v + v02 + v*v02/v01)
		b6 = -(xbar - xbar01)/(v + v01 + v*v01/v03) - (xbar - xbar03)/(v + v03 + v*v03/v01)
		b7 = -(xbar - xbar02)/(v + v02 + v*v02/v03) - (xbar - xbar03)/(v + v03 + v*v03/v02)
		b8 = -(xbar - xbar01)/(v + v01 + v*v01*(1/v02 + 1/v03)) - (xbar - xbar02)/(v + v02 + v*v02*(1/v01 + 1/v03)) - (xbar - xbar03)/(v + v03 + v*v03*(1/v01 + 1/v02))

		#Derivative of weight portion for delta method variance calculation (i.e., c1=the derivative of a1 wrt xbar)
		c1 = (w1*(b1-b1)+w2*(b2-b1)+w3*(b3-b1)+w4*(b4-b1)+w5*(b5-b1)+w6*(b6-b1)+w7*(b7-b1)+w8*(b8-b1))/w1
		c2 = (w1*(b1-b2)+w2*(b2-b2)+w3*(b3-b2)+w4*(b4-b2)+w5*(b5-b2)+w6*(b6-b2)+w7*(b7-b2)+w8*(b8-b2))/w2
		c3 = (w1*(b1-b3)+w2*(b2-b3)+w3*(b3-b3)+w4*(b4-b3)+w5*(b5-b3)+w6*(b6-b3)+w7*(b7-b3)+w8*(b8-b3))/w3
		c4 = (w1*(b1-b4)+w2*(b2-b4)+w3*(b3-b4)+w4*(b4-b4)+w5*(b5-b4)+w6*(b6-b4)+w7*(b7-b4)+w8*(b8-b4))/w4
		c5 = (w1*(b1-b5)+w2*(b2-b5)+w3*(b3-b5)+w4*(b4-b5)+w5*(b5-b5)+w6*(b6-b5)+w7*(b7-b5)+w8*(b8-b5))/w5
		c6 = (w1*(b1-b6)+w2*(b2-b6)+w3*(b3-b6)+w4*(b4-b6)+w5*(b5-b6)+w6*(b6-b6)+w7*(b7-b6)+w8*(b8-b6))/w6
		c7 = (w1*(b1-b7)+w2*(b2-b7)+w3*(b3-b7)+w4*(b4-b7)+w5*(b5-b7)+w6*(b6-b7)+w7*(b7-b7)+w8*(b8-b7))/w7
		c8 = (w1*(b1-b8)+w2*(b2-b8)+w3*(b3-b8)+w4*(b4-b8)+w5*(b5-b8)+w6*(b6-b8)+w7*(b7-b8)+w8*(b8-b8))/w8

		#Derivative of the posterior mean for delta method variance calculation
		dM1 = 1
		dM2 = v01/(v+v01)
		dM3 = v02/(v+v02)
		dM4 = v03/(v+v03)
		dM5 = v01*v02/(v01*v02 + v*v02 + v*v01)
		dM6 = v01*v03/(v01*v03 + v*v03 + v*v01)
		dM7 = v02*v03/(v02*v03 + v*v03 + v*v02)
		dM8 = v01*v02*v03/(v01*v02*v03 + v*v02*v03 + v*v01*v03 + v*v01*v02)

		#Derivative of g for delta method variance calculation
		g1 = (a1*dM1 - M1*c1) / a1^2
		g2 = (a2*dM2 - M2*c2) / a2^2
		g3 = (a3*dM3 - M3*c3) / a3^2
		g4 = (a4*dM4 - M4*c4) / a4^2
		g5 = (a5*dM5 - M5*c5) / a5^2
		g6 = (a6*dM6 - M6*c6) / a6^2
		g7 = (a7*dM7 - M7*c7) / a7^2
		g8 = (a8*dM8 - M8*c8) / a8^2

		gvec = c(g1,g2,g3,g4,g5,g6,g7,g8)
		gmat = gvec%*%t(gvec)

		mu.est <- w1*M1+w2*M2+w3*M3+w4*M4+w5*M5+w6*M6+w7*M7+w8*M8

		var.comp = rep(1,8)%*%gmat%*%rep(1,8) * v #accounting for covariance, variance of models
		bias.orig = ((w1*M1+w2*M2+w3*M3+w4*M4+w5*M5+w6*M6+w7*M7+w8*M8)-mu)
		mse = var.comp + bias.orig^2

		esss.orig = e1*w1+e2*w2+e3*w3+e4*w4+e5*w5+e6*w6+e7*w7+e8*w8

		post = t(rmultinom(1,30000,w))

		post.mod1 = rnorm(post[1,1],M1,sd=sqrt(V1))
		post.mod2 = rnorm(post[1,2],M2,sd=sqrt(V2))
		post.mod3 = rnorm(post[1,3],M3,sd=sqrt(V3))
		post.mod4 = rnorm(post[1,4],M4,sd=sqrt(V4))
		post.mod5 = rnorm(post[1,5],M5,sd=sqrt(V5))
		post.mod6 = rnorm(post[1,6],M6,sd=sqrt(V6))
		post.mod7 = rnorm(post[1,7],M7,sd=sqrt(V7))
		post.mod8 = rnorm(post[1,8],M8,sd=sqrt(V8))
		post.mu = c(post.mod1,post.mod2,post.mod3,post.mod4,post.mod5,post.mod6,post.mod7,post.mod8)
		hpd.bma_rep = boa.hpd(x=post.mu,alpha=.05)

		include.mu = (hpd.bma_rep[1] < mu) & (mu < hpd.bma_rep[2])
		diff.hpd = hpd.bma_rep[2] - hpd.bma_rep[1]

		est.mu <- c(est.mu, mu.est) #store simulation results for mean
		est.bias <- c(est.bias, bias.orig) #store simulation results for bias
		est.esss <- c(est.esss, esss.orig) #store simulation results for esss
		cov.ind = c(cov.ind, include.mu) #store simulation results for coverage
		hpd.width = c(hpd.width, diff.hpd) #store simulation results for HPD interval
	}

	var.part = var(est.mu) #estimate variance of estimated mu-hat from simulation
	mean.part = mean(est.mu) #estimate mu-hat
	bias.part = mean.part - mu #bias
	mse.part = bias.part^2 + var.part #MSE
	esss.part = median(est.esss) #ESSS
	esss.alt = NA #legacy variable	
	abserr.part = abs(bias.part) #absolute error of bias for summary stat
	cov.part = mean(cov.ind) #coverage
	hpd.part = mean(hpd.width) #HPD width

	mat = cbind(mu,mean.part,var.part,bias.part,mse.part,esss.part,esss.alt,abserr.part,cov.part,hpd.part)

	write.table(mat, file = paste0('MEM.estimators_',scen,'_',prior,'.txt'), append=TRUE, col.names=FALSE, row.names=FALSE)
#	return(mat)
}


calc.weights_MEM = function(xvec,svec,nvec,prior){
###function to calculate model weights for MEM approach with "correct" calculations which don't assume conditional independence
#xvec: means for sources
#svec: standard deviation for sources
#nvec: sample size for sources
#prior: prior to use for calculations

if(length(xvec)==4){
	xbar = xvec[1] #CENIC G
	xbar01 = xvec[2]#-0.15 #CENIC B
	xbar02 = xvec[3]#-4.24 #Quest Jr.
	xbar03 = xvec[4]#-7.08 #Quest Sr.

	s = svec[1]^2#6.79^2 #variance from CENIC G
	s01 = svec[2]^2#6.71^2
	s02 = svec[3]^2#9.02^2
	s03 = svec[4]^2#7.02^2

	n = nvec[1]#109
	n01 = nvec[2]#116
	n02 = nvec[3]#55
	n03 = nvec[4]#32

	v = s/n
	v01 = s01/n01
	v02 = s02/n02
	v03 = s03/n03

	###Calculating the weights
	#Writing out the marginal models
	m1 = sqrt(2*pi)^4 / sqrt(1/(v*v01*v02*v03))

	m2 = sqrt(2*pi)^3 / sqrt((1/v + 1/v01)*(1)/(v02*v03)) * exp(-0.5 * ((xbar-xbar01)^2/(v + v01)) )

	m3 = sqrt(2*pi)^3 / sqrt((1/v + 1/v02)*(1)/(v01*v03)) * exp(-0.5 * ((xbar-xbar02)^2/(v + v02)) )

	m4 = sqrt(2*pi)^3 / sqrt((1/v + 1/v03)*(1)/(v01*v02)) * exp(-0.5 * ((xbar-xbar03)^2/(v + v03)) )

	m5 = sqrt(2*pi)^2 / sqrt((1/v + 1/v01 + 1/v02)*(1/v03)) * exp(-0.5 * ( ((xbar-xbar01)^2/(v + v01 + (v*v01)/v02)) + ((xbar-xbar02)^2/(v + v02 + (v*v02)/v01)) + ((xbar01-xbar02)^2/(v01 + v02 + (v01*v02)/v)) ) )

	m6 = sqrt(2*pi)^2 / sqrt((1/v + 1/v01 + 1/v03)*(1/v02)) * exp(-0.5 * ( ((xbar-xbar01)^2/(v + v01 + (v*v01)/v03)) + ((xbar-xbar03)^2/(v + v03 + (v*v03)/v01)) + ((xbar01-xbar03)^2/(v01 + v03 + (v01*v03)/v)) ) )

	m7 = sqrt(2*pi)^2 / sqrt((1/v + 1/v02 + 1/v03)*(1/v01)) * exp(-0.5 * ( ((xbar-xbar02)^2/(v + v02 + (v*v02)/v03)) + ((xbar-xbar03)^2/(v + v03 + (v*v03)/v02)) + ((xbar02-xbar03)^2/(v02 + v03 + (v02*v03)/v)) ) )

	m8 = sqrt(2*pi) / sqrt(1/v + 1/v01 + 1/v02 + 1/v03) * exp(-0.5 * ( ((xbar-xbar01)^2/(v + v01 + (v*v01)*(v02^(-1)+v03^(-1)) )) + ((xbar-xbar02)^2/(v + v02 + (v*v02)*(v01^(-1)+v03^(-1)) )) + ((xbar-xbar03)^2/(v + v03 + (v*v03)*(v01^(-1)+v02^(-1)) )) + 
	((xbar01-xbar02)^2/(v01 + v02 + (v01*v02)*(v^(-1)+v03^(-1)) )) + ((xbar01-xbar03)^2/(v01 + v03 + (v01*v03)*(v^(-1)+v02^(-1)) )) + ((xbar02-xbar03)^2/(v02 + v03 + (v02*v03)*(v^(-1)+v01^(-1)) ))   ))

	if(prior=='pi_e'){
		prior1.1=prior2.1=prior3.1=prior1.0=prior2.0=prior3.0<-.5

		#note: don't need to use prior here since they're all equal to 1/2 and cancel out
		w1 = m1; w2 = m2; w3 = m3; w4 = m4; w5 = m5; w6 = m6; w7 = m7; w8 = m8
		m.sum.prior = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8
		q1 = m1 / m.sum.prior; q2 = m2 / m.sum.prior; q3 = m3 / m.sum.prior; q4 = m4 / m.sum.prior; q5 = m5 / m.sum.prior; q6 = m6 / m.sum.prior; q7 = m7 / m.sum.prior; q8 = m8 / m.sum.prior #weights for noninformative prior

	}else if(prior=='pi_n'){
		m1inv.noCurN = sqrt(1/(s*v01*v02*v03))/(2*pi)^(4/2) #only current study excluded
		m2inv.noCurN = sqrt((1/s + 1/v01)*(1/(v02*v03)))/(2*pi)^(3/2)
		m3inv.noCurN = sqrt((1/s + 1/v02)*(1/(v01*v03)))/(2*pi)^(3/2)
		m4inv.noCurN = sqrt((1/s + 1/v03)*(1/(v01*v02)))/(2*pi)^(3/2)
		m5inv.noCurN = sqrt((1/s + 1/v01 + 1/v02)*(1/v03))/(2*pi)^(2/2)
		m6inv.noCurN = sqrt((1/s + 1/v01 + 1/v03)*(1/v02))/(2*pi)^(2/2)
		m7inv.noCurN = sqrt((1/s + 1/v02 + 1/v03)*(1/v01))/(2*pi)^(2/2)
		m8inv.noCurN = sqrt(1/s + 1/v01 + 1/v02 + 1/v03)/(2*pi)^(1/2)

		prior1.1=m2inv.noCurN+m5inv.noCurN+m6inv.noCurN+m8inv.noCurN
		prior2.1=m3inv.noCurN+m5inv.noCurN+m7inv.noCurN+m8inv.noCurN
		prior3.1=m4inv.noCurN+m6inv.noCurN+m7inv.noCurN+m8inv.noCurN
		prior1.0=m1inv.noCurN+m3inv.noCurN+m4inv.noCurN+m7inv.noCurN
		prior2.0=m1inv.noCurN+m2inv.noCurN+m4inv.noCurN+m6inv.noCurN
		prior3.0=m1inv.noCurN+m2inv.noCurN+m3inv.noCurN+m5inv.noCurN

		m1pr.noCurN = prior1.0*prior2.0*prior3.0; m2pr.noCurN = prior1.1*prior2.0*prior3.0; m3pr.noCurN = prior1.0*prior2.1*prior3.0; m4pr.noCurN = prior1.0*prior2.0*prior3.1; m5pr.noCurN = prior1.1*prior2.1*prior3.0; m6pr.noCurN = prior1.1*prior2.0*prior3.1; m7pr.noCurN = prior1.0*prior2.1*prior3.1; m8pr.noCurN = prior1.1*prior2.1*prior3.1

		#Weight components:
		w1 = m1*m1pr.noCurN; w2 = m2*m2pr.noCurN; w3 = m3*m3pr.noCurN; w4 = m4*m4pr.noCurN; w5 = m5*m5pr.noCurN; w6 = m6*m6pr.noCurN; w7 = m7*m7pr.noCurN; w8 = m8*m8pr.noCurN #weight components for prior with no n
		m.sum.prior = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8

		#Final weights:
		q1 = w1 / m.sum.prior; q2 = w2 / m.sum.prior; q3 = w3 / m.sum.prior; q4 = w4 / m.sum.prior; q5 = w5 / m.sum.prior; q6 = w6 / m.sum.prior; q7 = w7 / m.sum.prior; q8 = w8 / m.sum.prior #weights for prior with no n

	}else{print('Only pi_e and pi_n currently supported.')}

	return(c(q1,q2,q3,q4,q5,q6,q7,q8))
}

if(length(xvec)==2){
	xbar = xvec[1] #CENIC G
	xbar01 = xvec[2]#-0.15 #CENIC B

	s = svec[1]^2#6.79^2 #variance from CENIC G
	s01 = svec[2]^2#6.71^2

	n = nvec[1]#109
	n01 = nvec[2]#116

	v = s/n
	v01 = s01/n01

	###Calculating the weights
	#Writing out the marginal models
	m1 = sqrt(2*pi)^2 / sqrt(1/(v*v01))

	m2 = sqrt(2*pi)^1 / sqrt((1/v + 1/v01)) * exp(-0.5 * ((xbar-xbar01)^2/(v + v01)) )

	if(prior=='pi_e'){
		w1 <- 0.5*m1; w2 <- 0.5*m2
		m.sum.prior <- w1 + w2
		
		q1 <- w1/m.sum.prior; q2 <- w2/m.sum.prior
	}else if(prior=='pi_n'){
		m1inv.noCurN = sqrt(1/(s*v01))/(2*pi)^(2/2) #only current study excluded
		m2inv.noCurN = sqrt(1/s + 1/v01)/(2*pi)^(1/2)

		prior1.1 = m2inv.noCurN
		prior1.0 = m1inv.noCurN

		m1pr <- prior1.0; m2pr <- prior1.1

		w1 = m1*m1pr; w2 = m2*m2pr 
		m.sum.prior = w1 + w2

		#Final weights:
		q1 = w1 / m.sum.prior; q2 = w2 / m.sum.prior 
	}else{print('Only pi_e and pi_n currently supported.')}

	return(c(q1,q2))
}

}

#############################################
### MEM functions related to convergence of models
#############################################

weight.tab.BMA <- function(xvec,svec,nvec,ncur,prior,subtext=NULL,plot=F,cur.n_inf=T){
###Function to create table showing model weights and plot trajectories for BMA approach
#plot=T indicates a plot of table should also be produced
#cur.n_inf indicates if all sample sizes should go to infinity (F) or just the current study (T)
#subtext is a text string to optionally place as subtitle on plot

	if(cur.n_inf==T){
		tab <- sapply(1:length(ncur), function(i) calc.weights_BMA(x=xvec,S=svec,N=c(ncur[i],nvec),prior=prior))
		xlab.asd <- 'Current Study Sample Size'
	}else{
		tab <- sapply(1:length(ncur), function(i) calc.weights_BMA(x=xvec,S=svec,N=rep(ncur[i],4),prior=prior))
		xlab.asd <- 'Sample Size for all Studies'
	}
	colnames(tab) <- ncur
	rownames(tab) <- c('M1','M2','M3','M4','M5','M6','M7','M8')
	
	if(plot==T){
		cols = rainbow(7)
		plot(x=ncur,y=tab[1,],type='l',ylim=c(0,1),ylab='Posterior Model Weight',xlab=xlab.asd);text(x=ncur[length(ncur)-50],y=tab[1,length(ncur)-1],labels='M1')
		mtext(subtext,line=.5)
		lines(x=ncur,y=tab[2,],col=cols[1],lty=2);text(x=ncur[length(ncur)-100],y=tab[2,length(ncur)-100],labels='M2',col=cols[1])
		lines(x=ncur,y=tab[3,],col=cols[2],lty=2);text(x=ncur[length(ncur)-150],y=tab[3,length(ncur)-150],labels='M3',col=cols[2])
		lines(x=ncur,y=tab[4,],col=cols[3],lty=2);text(x=ncur[length(ncur)-200],y=tab[4,length(ncur)-200],labels='M4',col=cols[3])
		lines(x=ncur,y=tab[5,],col=cols[4],lty=6);text(x=ncur[length(ncur)-250],y=tab[5,length(ncur)-250],labels='M5',col=cols[4])
		lines(x=ncur,y=tab[6,],col=cols[5],lty=4);text(x=ncur[length(ncur)-300],y=tab[6,length(ncur)-300],labels='M6',col=cols[5])
		lines(x=ncur,y=tab[7,],col=cols[6],lty=4);text(x=ncur[length(ncur)-350],y=tab[7,length(ncur)-350],labels='M7',col=cols[6])
		lines(x=ncur,y=tab[8,],col=cols[7],lty=5);text(x=ncur[length(ncur)-400],y=tab[8,length(ncur)-400],labels='M8',col=cols[7])
	}
	tab.return <- tab[,which(colnames(tab)%in%c(1,10,100,1000,10000,100000,1000000,10000000))]
	return(tab.return)
}


weight.tab.MEM <- function(xvec,svec,nvec,ncur,prior,subtext=NULL,plot=F,cur.n_inf=T){
###Function to create table showing model weights and plot trajectories for MEM approach
#plot=T indicates a plot of table should also be produced
#cur.n_inf=T indicates that only the current sample size goes to infinity, F means all studies
#subtext is a text string to optionally place as subtitle on plot

	if(cur.n_inf==T){
		tab <- sapply(1:length(ncur), function(i) calc.weights_MEM(xvec=xvec,svec=svec,nvec=c(ncur[i],nvec),prior=prior))
		xlab.asd <- 'Current Study Sample Size'
	}else{
		tab <- sapply(1:length(ncur), function(i) calc.weights_MEM(xvec=xvec,svec=svec,nvec=rep(ncur[i],4),prior=prior))
		xlab.asd <- 'Sample Size for all Studies'
	}
	colnames(tab) <- ncur
	rownames(tab) <- c('M1','M2','M3','M4','M5','M6','M7','M8')
	
	if(plot==T){
		cols = rainbow(7)
		plot(x=ncur,y=tab[1,],type='l',ylim=c(0,1),ylab='Posterior Model Weight',xlab=xlab.asd);text(x=ncur[length(ncur)-50],y=tab[1,length(ncur)-1],labels='M1')
		#main='drBMA'
		mtext(subtext,line=.5)
		lines(x=ncur,y=tab[2,],col=cols[1],lty=2);text(x=ncur[length(ncur)-100],y=tab[2,length(ncur)-100],labels='M2',col=cols[1])
		lines(x=ncur,y=tab[3,],col=cols[2],lty=2);text(x=ncur[length(ncur)-150],y=tab[3,length(ncur)-150],labels='M3',col=cols[2])
		lines(x=ncur,y=tab[4,],col=cols[3],lty=2);text(x=ncur[length(ncur)-200],y=tab[4,length(ncur)-200],labels='M4',col=cols[3])
		lines(x=ncur,y=tab[5,],col=cols[4],lty=6);text(x=ncur[length(ncur)-250],y=tab[5,length(ncur)-250],labels='M5',col=cols[4])
		lines(x=ncur,y=tab[6,],col=cols[5],lty=4);text(x=ncur[length(ncur)-300],y=tab[6,length(ncur)-300],labels='M6',col=cols[5])
		lines(x=ncur,y=tab[7,],col=cols[6],lty=4);text(x=ncur[length(ncur)-350],y=tab[7,length(ncur)-350],labels='M7',col=cols[6])
		lines(x=ncur,y=tab[8,],col=cols[7],lty=5);text(x=ncur[length(ncur)-400],y=tab[8,length(ncur)-400],labels='M8',col=cols[7])
	}
	tab.return <- tab[,which(colnames(tab)%in%c(1,10,100,1000,10000,100000,1000000,10000000))]
	return(tab.return)
}

calc.weights_BMA = function(prior,x,S,N){
###Function to calculate BMA weights for case with 1 or 3 supplemental sources
#x: means for primary and then supplemental sources
#S: standard deviation for primary and then supplemental sources
#N: sample size for primary and then supplemental sources
#prior: chose which prior to fit model on ('pi_e','pi_n','minv.N','minv.noN')

if(length(S)==4){
	xbar = x[1] 
	xbar01 = x[2] 
	xbar02 = x[3] 
	xbar03 = x[4] 

	s = S[1]^2 
	s01 = S[2]^2
	s02 = S[3]^2
	s03 = S[4]^2

	n = N[1]
	n01 = N[2]
	n02 = N[3]
	n03 = N[4]

	v = s/n
	v01 = s01/n01
	v02 = s02/n02
	v03 = s03/n03

	###Calculating the weights
	#Writing out the marginal models
	m1 = sqrt(2*pi)^4 / sqrt(1/(v*v01*v02*v03))

	m2 = sqrt(2*pi)^3 / sqrt((1/v + 1/v01)*(1)/(v02*v03)) * exp(-0.5 * ((xbar-xbar01)^2/(v + v01)) )

	m3 = sqrt(2*pi)^3 / sqrt((1/v + 1/v02)*(1)/(v01*v03)) * exp(-0.5 * ((xbar-xbar02)^2/(v + v02)) )

	m4 = sqrt(2*pi)^3 / sqrt((1/v + 1/v03)*(1)/(v01*v02)) * exp(-0.5 * ((xbar-xbar03)^2/(v+ v03)) )

	m5 = sqrt(2*pi)^2 / sqrt((1/v + 1/v01 + 1/v02)*(1/v03)) * exp(-0.5 * ( ((xbar-xbar01)^2/(v + v01 + (v*v01)/v02)) + ((xbar-xbar02)^2/(v + v02 + (v*v02)/v01)) + ((xbar01-xbar02)^2/(v01 + v02 + (v01*v02)/v)) ) )

	m6 = sqrt(2*pi)^2 / sqrt((1/v + 1/v01 + 1/v03)*(1/v02)) * exp(-0.5 * ( ((xbar-xbar01)^2/(v + v01 + (v*v01)/v03)) + ((xbar-xbar03)^2/(v + v03 + (v*v03)/v01)) + ((xbar01-xbar03)^2/(v01 + v03 + (v01*v03)/v)) ) )

	m7 = sqrt(2*pi)^2 / sqrt((1/v + 1/v02 + 1/v03)*(1/v01)) * exp(-0.5 * ( ((xbar-xbar02)^2/(v + v02 + (v*v02)/v03)) + ((xbar-xbar03)^2/(v + v03 + (v*v03)/v02)) + ((xbar02-xbar03)^2/(v02 + v03 + (v02*v03)/v)) ) )

	m8 = sqrt(2*pi) / sqrt(1/v + 1/v01 + 1/v02 + 1/v03) * exp(-0.5 * ( ((xbar-xbar01)^2/(v + v01 + (v*v01)*(v02^(-1)+v03^(-1)) )) + ((xbar-xbar02)^2/(v + v02 + (v*v02)*(v01^(-1)+v03^(-1)) )) + ((xbar-xbar03)^2/(v + v03 + (v*v03)*(v01^(-1)+v02^(-1)) )) + 
	((xbar01-xbar02)^2/(v01 + v02 + (v01*v02)*(v^(-1)+v03^(-1)) )) + ((xbar01-xbar03)^2/(v01 + v03 + (v01*v03)*(v^(-1)+v02^(-1)) )) + ((xbar02-xbar03)^2/(v02 + v03 + (v02*v03)*(v^(-1)+v01^(-1)) ))   ))

	m.sum = m1+m2+m3+m4+m5+m6+m7+m8

	#Priors
	if(prior=='pi_e'){
		w1 = m1*(1/8); w2 = m2*(1/8); w3 = m3*(1/8); w4 = m4*(1/8); w5 = m5*(1/8); w6 = m6*(1/8); w7 = m7*(1/8); w8 = m8*(1/8) #uninform prior
		m.sum.prior = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8
		q1 = w1 / m.sum.prior; q2 = w2 / m.sum.prior; q3 = w3 / m.sum.prior; q4 = w4 / m.sum.prior; q5 = w5 / m.sum.prior; q6 = w6 / m.sum.prior; q7 = w7 / m.sum.prior; q8 = w8 / m.sum.prior #weights for noninformative prior
	}else if(prior=='pi_n'){
		m1inv.noCurN = sqrt(1/(s*v01*v02*v03))/(2*pi)^(4/2) #only current study N excluded
		m2inv.noCurN = sqrt((1/s + 1/v01)*(1/(v02*v03)))/(2*pi)^(3/2)
		m3inv.noCurN = sqrt((1/s + 1/v02)*(1/(v01*v03)))/(2*pi)^(3/2)
		m4inv.noCurN = sqrt((1/s + 1/v03)*(1/(v01*v02)))/(2*pi)^(3/2)
		m5inv.noCurN = sqrt((1/s + 1/v01 + 1/v02)*(1/v03))/(2*pi)^(2/2)
		m6inv.noCurN = sqrt((1/s + 1/v01 + 1/v03)*(1/v02))/(2*pi)^(2/2)
		m7inv.noCurN = sqrt((1/s + 1/v02 + 1/v03)*(1/v01))/(2*pi)^(2/2)
		m8inv.noCurN = sqrt(1/s + 1/v01 + 1/v02 + 1/v03)/(2*pi)^(1/2)

		m1pr.noCurN = m1inv.noCurN; m2pr.noCurN = m2inv.noCurN; m3pr.noCurN = m3inv.noCurN; m4pr.noCurN = m4inv.noCurN; m5pr.noCurN = m5inv.noCurN; m6pr.noCurN = m6inv.noCurN; m7pr.noCurN = m7inv.noCurN; m8pr.noCurN = m8inv.noCurN

		#Weight components:
		w1 = m1*m1pr.noCurN; w2 = m2*m2pr.noCurN; w3 = m3*m3pr.noCurN; w4 = m4*m4pr.noCurN; w5 = m5*m5pr.noCurN; w6 = m6*m6pr.noCurN; w7 = m7*m7pr.noCurN; w8 = m8*m8pr.noCurN #weight components for prior with no n
		m.sum.prior = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8

		#Final weights:
		q1 = m1*m1pr.noCurN / m.sum.prior; q2 = m2*m2pr.noCurN / m.sum.prior; q3 = m3*m3pr.noCurN / m.sum.prior; q4 = m4*m4pr.noCurN / m.sum.prior; q5 = m5*m5pr.noCurN / m.sum.prior; q6 = m6*m6pr.noCurN / m.sum.prior; q7 = m7*m7pr.noCurN / m.sum.prior; q8 = m8*m8pr.noCurN / m.sum.prior #weights for prior with no n
	}else{print('Prior not recognized, unable to calculate model weights.')}

	return(c(q1,q2,q3,q4,q5,q6,q7,q8))
}
if(length(S)==2){
	xbar = x[1] 
	xbar01 = x[2] 

	s = S[1]^2 
	s01 = S[2]^2

	n = N[1]
	n01 = N[2]

	v = s/n
	v01 = s01/n01

	###Calculating the weights
	#Writing out the marginal models
	m1 = sqrt(2*pi)^2 / sqrt(1/(v*v01))

	m2 = sqrt(2*pi)^1 / sqrt(1/v + 1/v01) * exp(-0.5 * ((xbar-xbar01)^2/(v + v01)) )

	m.sum = m1+m2

	#Priors
	if(prior=='pi_e'){
		w1 = m1*.5; w2 = m2*.5 #uniform prior
		m.sum.prior = w1 + w2 
		q1 = m1 / m.sum.prior; q2 = m2 / m.sum.prior #weights for noninformative prior
	}else if(prior=='pi_n'){
		m1inv.noCurN = sqrt(1/(s*v01))/(2*pi)^(2/2) #only current study N excluded
		m2inv.noCurN = sqrt(1/s + 1/v01)/(2*pi)^(1/2)

		m1pr.noCurN = m1inv.noCurN; m2pr.noCurN = m2inv.noCurN

		#Weight components:
		w1 = m1*m1pr.noCurN; w2 = m2*m2pr.noCurN #weight components for prior with no n
		m.sum.prior = w1 + w2

		#Final weights:
		q1 = m1*m1pr.noCurN / m.sum.prior; q2 = m2*m2pr.noCurN / m.sum.prior #weights for prior with no n
	}else{print('Prior not recognized, unable to calculate model weights.')}

	return(c(q1,q2))
}
}

#############################################
### Functions related to CENIC results
#############################################

summary_calc = function(x,S,N,prior,type){
###x: vector mean values for studies
###S: vector of standard deviations for studies
###N: sample size for studies
###prior: prior to use for calculation
###type: 'BMA' for BMA framework, 'MEM' for MEM framework

if(length(S)==4){
	xbar = x[1]
	xbar01 = x[2] #CENIC B
	xbar02 = x[3] #Quest Jr.
	xbar03 = x[4] #Quest Sr.

	s = S[1]^2 #variance
	s01 = S[2]^2
	s02 = S[3]^2
	s03 = S[4]^2

	n = N[1] #sample sizes
	n01 = N[2]
	n02 = N[3]
	n03 = N[4]

	v = s/n
	v01 = s01/n01
	v02 = s02/n02
	v03 = s03/n03
	
	###Posterior components
	##Means for each model
	M1 = xbar
	M2 = (v01*xbar + v*xbar01)/(v+v01)
	M3 = (v02*xbar + v*xbar02)/(v+v02)
	M4 = (v03*xbar + v*xbar03)/(v+v03)
	M5 = (v01*v02*xbar + v*v02*xbar01 + v*v01*xbar02)/(v01*v02 + v*v02 + v*v01)
	M6 = (v01*v03*xbar + v*v03*xbar01 + v*v01*xbar03)/(v01*v03 + v*v03 + v*v01)
	M7 = (v02*v03*xbar + v*v03*xbar02 + v*v02*xbar03)/(v02*v03 + v*v03 + v*v02)
	M8 = (v01*v02*v03*xbar + v*v02*v03*xbar01 + v*v01*v03*xbar02 + v*v01*v02*xbar03)/(v01*v02*v03 + v*v02*v03 + v*v01*v03 + v*v01*v02)

	##Variances for each model
	V1 = v
	V2 = (1/v + 1/v01)^-1
	V3 = (1/v + 1/v02)^-1
	V4 = (1/v + 1/v03)^-1
	V5 = (1/v + 1/v01 + 1/v02)^-1
	V6 = (1/v + 1/v01 + 1/v03)^-1
	V7 = (1/v + 1/v02 + 1/v03)^-1
	V8 = (1/v + 1/v01 + 1/v02 + 1/v03)^-1

	##Effective historical sample size for each model
	e1 = (s/V1) - n
	e2 = (s/V2) - n
	e3 = (s/V3) - n
	e4 = (s/V4) - n
	e5 = (s/V5) - n
	e6 = (s/V6) - n
	e7 = (s/V7) - n
	e8 = (s/V8) - n

	if(type=='MEM'){
		w = calc.weights_MEM(xvec=c(xbar,x[2],x[3],x[4]),svec=c(S[1],S[2],S[3],S[4]),nvec=c(N[1],N[2],N[3],N[4]),prior=prior)
	}else if(type=='BMA'){
		w = calc.weights.multinomial(x=c(xbar,x[2],x[3],x[4]),S=c(S[1],S[2],S[3],S[4]),N=c(N[1],N[2],N[3],N[4]),prior=prior)
	}else{print('Prior framework not supported, use MEM or BMA.')}

	w[w < 1.630888e-293] <-  1.630888e-293
	w1<-w[1]; w2<-w[2]; w3<-w[3]; w4<-w[4]; w5<-w[5]; w6<-w[6]; w7<-w[7]; w8<-w[8]

	###Calculate MSE
	##Variance component
	#Weight denominator needed for delta method g function calculation (weight for each model):
	##Weight is w1 = 1/a1, calculations as a1 for ease of inclusion with quotient rule for derivatives
	a1 = (w1+w2+w3+w4+w5+w6+w7+w8)/w1 #note that 1/a1 equals the weight for model 1, and 1/a2 for model 2, etc.
	a2 = (w1+w2+w3+w4+w5+w6+w7+w8)/w2
	a3 = (w1+w2+w3+w4+w5+w6+w7+w8)/w3
	a4 = (w1+w2+w3+w4+w5+w6+w7+w8)/w4
	a5 = (w1+w2+w3+w4+w5+w6+w7+w8)/w5
	a6 = (w1+w2+w3+w4+w5+w6+w7+w8)/w6
	a7 = (w1+w2+w3+w4+w5+w6+w7+w8)/w7
	a8 = (w1+w2+w3+w4+w5+w6+w7+w8)/w8

	#Derivative of exponential part of each w_i weight component for delta method variance calculation
	##Note that the exp() part is left off because it is contained within the w_i part included in the next step
	##(i.e., (w1+w2+w3+...+w8)/w2 will contain all the necessary info besides the deriv of the exp part calculated for b here)
	b1 = 0
	b2 = -(xbar - xbar01)/(v + v01)
	b3 = -(xbar - xbar02)/(v + v02)
	b4 = -(xbar - xbar03)/(v + v03)
	b5 = -(xbar - xbar01)/(v + v01 + v*v01/v02) - (xbar - xbar02)/(v + v02 + v*v02/v01)
	b6 = -(xbar - xbar01)/(v + v01 + v*v01/v03) - (xbar - xbar03)/(v + v03 + v*v03/v01)
	b7 = -(xbar - xbar02)/(v + v02 + v*v02/v03) - (xbar - xbar03)/(v + v03 + v*v03/v02)
	b8 = -(xbar - xbar01)/(v + v01 + v*v01*(1/v02 + 1/v03)) - (xbar - xbar02)/(v + v02 + v*v02*(1/v01 + 1/v03)) - (xbar - xbar03)/(v + v03 + v*v03*(1/v01 + 1/v02))

	#Derivative of weight portion for delta method variance calculation (i.e., c1=the derivative of a1 wrt xbar)
	c1 = (w1*(b1-b1)+w2*(b2-b1)+w3*(b3-b1)+w4*(b4-b1)+w5*(b5-b1)+w6*(b6-b1)+w7*(b7-b1)+w8*(b8-b1))/w1
	c2 = (w1*(b1-b2)+w2*(b2-b2)+w3*(b3-b2)+w4*(b4-b2)+w5*(b5-b2)+w6*(b6-b2)+w7*(b7-b2)+w8*(b8-b2))/w2
	c3 = (w1*(b1-b3)+w2*(b2-b3)+w3*(b3-b3)+w4*(b4-b3)+w5*(b5-b3)+w6*(b6-b3)+w7*(b7-b3)+w8*(b8-b3))/w3
	c4 = (w1*(b1-b4)+w2*(b2-b4)+w3*(b3-b4)+w4*(b4-b4)+w5*(b5-b4)+w6*(b6-b4)+w7*(b7-b4)+w8*(b8-b4))/w4
	c5 = (w1*(b1-b5)+w2*(b2-b5)+w3*(b3-b5)+w4*(b4-b5)+w5*(b5-b5)+w6*(b6-b5)+w7*(b7-b5)+w8*(b8-b5))/w5
	c6 = (w1*(b1-b6)+w2*(b2-b6)+w3*(b3-b6)+w4*(b4-b6)+w5*(b5-b6)+w6*(b6-b6)+w7*(b7-b6)+w8*(b8-b6))/w6
	c7 = (w1*(b1-b7)+w2*(b2-b7)+w3*(b3-b7)+w4*(b4-b7)+w5*(b5-b7)+w6*(b6-b7)+w7*(b7-b7)+w8*(b8-b7))/w7
	c8 = (w1*(b1-b8)+w2*(b2-b8)+w3*(b3-b8)+w4*(b4-b8)+w5*(b5-b8)+w6*(b6-b8)+w7*(b7-b8)+w8*(b8-b8))/w8

	#Derivative of the posterior mean for delta method variance calculation
	dM1 = 1
	dM2 = v01/(v+v01)
	dM3 = v02/(v+v02)
	dM4 = v03/(v+v03)
	dM5 = v01*v02/(v01*v02 + v*v02 + v*v01)
	dM6 = v01*v03/(v01*v03 + v*v03 + v*v01)
	dM7 = v02*v03/(v02*v03 + v*v03 + v*v02)
	dM8 = v01*v02*v03/(v01*v02*v03 + v*v02*v03 + v*v01*v03 + v*v01*v02)

	#Derivative of g for delta method variance calculation
	g1 = (a1*dM1 - M1*c1) / a1^2
	g2 = (a2*dM2 - M2*c2) / a2^2
	g3 = (a3*dM3 - M3*c3) / a3^2
	g4 = (a4*dM4 - M4*c4) / a4^2
	g5 = (a5*dM5 - M5*c5) / a5^2
	g6 = (a6*dM6 - M6*c6) / a6^2
	g7 = (a7*dM7 - M7*c7) / a7^2
	g8 = (a8*dM8 - M8*c8) / a8^2

	gvec = c(g1,g2,g3,g4,g5,g6,g7,g8)
	gmat = gvec%*%t(gvec)

	mean.est <- w1*M1+w2*M2+w3*M3+w4*M4+w5*M5+w6*M6+w7*M7+w8*M8
	var.est = rep(1,8)%*%gmat%*%rep(1,8) * v #accounting for covariance, variance of models
	esss.est <- e1*w1 + e2*w2 + e3*w3 + e4*w4 + e5*w5 + e6*w6 + e7*w7 + e8*w8

	est.vec = c(mean.est, var.est, esss.est,w)
}
if(length(S)==2){
	xbar = x[1]
	xbar01 = x[2] #CENIC B

	s = S[1]^2 #variance
	s01 = S[2]^2

	n = N[1] #sample sizes
	n01 = N[2]

	v = s/n
	v01 = s01/n01
	
	###Posterior components
	##Means for each model
	M1 = xbar
	M2 = (v01*xbar + v*xbar01)/(v+v01)

	##Variances for each model
	V1 = v
	V2 = (1/v + 1/v01)^-1

	##Effective historical sample size for each model
	e1 = (s/V1) - n
	e2 = (s/V2) - n

	if(type=='MEM'){
		w = calc.weights_MEM(xvec=c(xbar,x[2]),svec=c(S[1],S[2]),nvec=c(N[1],N[2]),prior=prior)
	}else if(type=='BMA'){
		w = calc.weights.multinomial(x=c(xbar,x[2]),S=c(S[1],S[2]),N=c(N[1],N[2]),prior=prior)
	}else{print('Prior framework not supported, use MEM or BMA.')}

	w[w < 1.630888e-293] <-  1.630888e-293
	w1<-w[1]; w2<-w[2]

	###Calculate MSE
	##Variance component
	#Weight denominator needed for delta method g function calculation (weight for each model):
	##Weight is w1 = 1/a1, calculations as a1 for ease of inclusion with quotient rule for derivatives
	a1 = (w1+w2)/w1 #note that 1/a1 equals the weight for model 1, and 1/a2 for model 2, etc.
	a2 = (w1+w2)/w2

	#Derivative of exponential part of each w_i weight component for delta method variance calculation
	##Note that the exp() part is left off because it is contained within the w_i part included in the next step
	##(i.e., (w1+w2+w3+...+w8)/w2 will contain all the necessary info besides the deriv of the exp part calculated for b here)
	b1 = 0
	b2 = -(xbar - xbar01)/(v + v01)

	#Derivative of weight portion for delta method variance calculation (i.e., c1=the derivative of a1 wrt xbar)
	c1 = (w1*(b1-b1)+w2*(b2-b1))/w1
	c2 = (w1*(b1-b2)+w2*(b2-b2))/w2

	#Derivative of the posterior mean for delta method variance calculation
	dM1 = 1
	dM2 = v01/(v+v01)

	#Derivative of g for delta method variance calculation
	g1 = (a1*dM1 - M1*c1) / a1^2
	g2 = (a2*dM2 - M2*c2) / a2^2

	gvec = c(g1,g2)
	gmat = gvec%*%t(gvec)

	mean.est <- w1*M1+w2*M2
	var.est <- rep(1,2)%*%gmat%*%rep(1,2) * v #accounting for covariance, variance of models
	esss.est <- e1*w1 + e2*w2

	est.vec = c(mean.est, var.est, esss.est,w)
}

	return(est.vec)
}






