library(xtable)

source('MEM_Functions.R') #set source code file

####We note that our code requires only the sample mean, sample standard deviation and sample size from each source####
####basic summary statistics can be found below###

#Means for CENIC-p1 arms and Quest studies
g.mean = -0.2274181
b.mean =  -0.1477073
c.mean = 5.896042
e.mean = 7.33166
jr.mean = -4.244408
sr.mean = -7.081956

#Standard deivations for CENIC-p1 arms and Quest studies
g.sd = 6.78897
b.sd = 6.705689
c.sd = 9.145066
e.sd = 8.382593
jr.sd = 9.022088
sr.sd = 7.020726

#N's (sample sizes) for CENIC-p1 arms and Quest studies
g.n <- 109
b.n <- 116
c.n <- 110
e.n <- 112
jr.n <- 55
sr.n <- 32

###z-test function
z.test = function(mean.ctrl, mean.trt, se.ctrl, se.trt){
	z <- abs(mean.ctrl - mean.trt)/sqrt(se.ctrl + se.trt)
	return( pnorm(z,lower.tail=F) )
}

#store results
res.mat <- matrix(nrow=15,ncol=5)
rownames(res.mat) <- c('w1.int','w2.int','w3.int','w4.int','w5.int','w6.int','w7.int','w8.int','w1.con','w2.con','trt','esss.trt','ctrl','esss.ctrl','trt.effect')
colnames(res.mat) <- c('No Borrowing','MEM: pi_e','MEM: pi_n','Comm Pr','Meta-Analysis')

#no borrowing:
res.mat[,1] <- c(c(1,rep(0,7),1,0),paste0(round(g.mean,3),' (',round(g.sd/sqrt(g.n),3),')'),0,paste0(round(c.mean,3),' (',round(c.sd/sqrt(c.n),3),')'),0, paste0(round(g.mean-c.mean,3),' (',round(sqrt(g.sd^2/g.n + c.sd^2/c.n),3),')'))

#MEM
priors <- c('pi_e','pi_n')
for(i in 1:length(priors)){
	pr.use <- priors[i]
	ctrl <- summary_calc(x=c(c.mean,e.mean),S=c(c.sd,e.sd),N=c(c.n,e.n),prior=pr.use,type='MEM')
	trt <- summary_calc(x=c(g.mean,b.mean,jr.mean,sr.mean),S=c(g.sd,b.sd,jr.sd,sr.sd),N=c(g.n,b.n,jr.n,sr.n),prior=pr.use,type='MEM')
	res.mat[,(i+1)] <- c(round(trt[4:11],3),round(ctrl[4:5],3),paste0(round(trt[1],3),' (',round(sqrt(trt[2]),3),')'),round(trt[3],1),paste0(round(ctrl[1],3),' (',round(sqrt(ctrl[2]),3),')'),round(ctrl[3],1), paste0(round(trt[1]-ctrl[1],3),' (',round(sqrt(trt[2]+ctrl[2]),3),')'))
}

#commensurate prior:
ctrl <- commensurateprior_calc(x=c(c.mean,e.mean),S=c(c.sd,e.sd),N=c(c.n,e.n))
trt <- commensurateprior_calc(x=c(g.mean,b.mean,jr.mean,sr.mean),S=c(g.sd,b.sd,jr.sd,sr.sd),N=c(g.n,b.n,jr.n,sr.n))
res.mat[,4] <- c(rep(NA,10),paste0(round(trt[1],3),' (',round(sqrt(trt[2]),3),')'),round(trt[3],1),paste0(round(ctrl[1],3),' (',round(sqrt(ctrl[2]),3),')'),round(ctrl[3],1), paste0(round(trt[1]-ctrl[1],3),' (',round(sqrt(trt[2]+ctrl[2]),3),')'))

#meta-analysis:
##setwd() #set directory where BUGS programs are stored to run standard hierarchical model calculations

ctrl <- shm_calc(x=c(c.mean,e.mean),S=c(c.sd,e.sd),N=c(c.n,e.n),tau.up=50)
trt <- shm_calc(x=c(g.mean,b.mean,jr.mean,sr.mean),S=c(g.sd,b.sd,jr.sd,sr.sd),N=c(g.n,b.n,jr.n,sr.n),tau.up=50)
res.mat[,5] <- c(rep(NA,10),paste0(round(trt[1],3),' (',round(sqrt(trt[2]),3),')'),round(trt[3],1),paste0(round(ctrl[1],3),' (',round(sqrt(ctrl[2]),3),')'),round(ctrl[3],1), paste0(round(trt[1]-ctrl[1],3),' (',round(sqrt(trt[2]+ctrl[2]),3),')'))

