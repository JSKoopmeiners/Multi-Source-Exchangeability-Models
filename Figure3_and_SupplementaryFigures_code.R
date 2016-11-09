
###set appropriate working directory for R###

setwd() 

####Code to create plots of various summaries of performance as a function of the true mu###
####The code is currently set up to create plots for the summaries considered in the manuscript###
####Data for simulations considered in the manuscript can be found in the Simulation Results sub-folder####
####simulation code that can be edited to accomodate new scenarios can be found in MEM_Simulations.R###

####Figure 3 from paper is ESSS as a function of mu - summary = 6###
####Supplementary figures include: bias = 4, mse = 5 and coverage = 9###
####Figures can also be created for var = 3, abs.error = 8 and hpd width = 10 but are not presented in the supplementary material or manuscript###

scen.plot <- 1   ####identify scenario that should be used for making plots: paper used scenarios 1 through 4###
summary <- 9   ####identify summary that should be considered, options include: bias=4, var=3, mse=5, esss=6, abs.error=8, coverage=9, hpd width=10###

prior.plot <- c('pi_n.txt','pi_e.txt')
MEM.list <- list()
for(i in 1:2){
	asd.MEM <- read.table(paste0('Simulation Results/MEM.estimators_scen',scen.plot,'_',prior.plot[i]),header=F)
	MEM.list[[i]] <- asd.MEM[order(asd.MEM[,1]),]
}

commpr <- read.table(paste0('Simulation Results/commprior.estimators_scen',scen.plot,'.txt'),header=F); commpr <- commpr[order(commpr[,1]),]
meta.app <- read.table(paste0('Simulation Results/shm.estimators_scen',scen.plot,'.txt'),header=F); meta.app[which(meta.app[,6]<0),6] <- 0; meta.app <- meta.app[order(meta.app[,1]),]


#plotting results
plot.title = c(paste0('Estimated Bias, Scenario ',scen.plot),paste0('Estimated Variance, Scenario ',scen.plot),paste0('Estimated MSE, Scenario ',scen.plot),paste0('Scenario ', scen.plot),paste0('Average Coverage, Scenario ', scen.plot))
plot.ylab = c('Bias','Variance','MSE','Median ESSS','Average Coverage')
pdf.lab = c('Bias','Variance','MSE','ESSS','Coverage')

sdvec = c(4,4,4,4)
Nvec = c(100,100,100,100)
if(scen.plot==1){
	xvec = c(-4,-4,-4)
}else if(scen.plot==2){
	xvec = c(-10,-10,2)
}else if(scen.plot==3){
	xvec = c(-10,-4,2)
}else if(scen.plot==4){
	xvec = c(-10,-9.25,2)
}else{print('No scenario given.')}

V1.679 = (sdvec[1]/sqrt(Nvec[1]))^2
xbar01 = xvec[1]
xbar02 = xvec[2]
xbar03 = xvec[3]
hline = c(0,V1.679,V1.679,NA,NA)

cvec = c(4,3,5,6,9) #column of results to use for plotting mean: bias=4, var=3, mse=5, esss=6, abs. error=8, coverage=9, hpd width=10
ylim.plot = list(c(-.3,.3), c(0,.45), c(0,.3), c(0,175), c(.8,1)) #y-axis bounds for plots

cols = c('#000000','#e41a1c','darkgreen','#4daf4a','blue','#ff7f00','orangered','#f781bf','gray65')


	i <- which(cvec == summary)

	c.use = cvec[i]
	plot(x=MEM.list[[1]][,1],y=MEM.list[[1]][,c.use],xlim=c(-15,6),col=cols[1],main=plot.title[i],ylab=plot.ylab[i],xlab=expression(mu),type='l',lwd=2,ylim=ylim.plot[[i]]) 
	abline(h=hline[i],col='lightgray')
	lines(x=MEM.list[[2]][,1],y=MEM.list[[2]][,c.use],lwd=2,col=cols[5],lty=6)
	lines(smooth.spline(x=meta.app[,1],y=meta.app[,c.use],spar=.5),lwd=2,lty=4,col=cols[7])
	lines(x=commpr[,1],y=commpr[,c.use],lwd=2,lty=5,col=cols[9])
	abline(v=c(xbar01,xbar02,xbar03),col='lightgray',lty=2)
	plot.legend <- c(expression(paste('MEM ',pi[n])),expression(paste('MEM ',pi[n^"l"])),expression(paste('MEM ',pi[e])),'CP','SHM')
	if((scen.plot == 1)){
	legend('topleft',col=cols[c(1,5,9,7)],lty=c(1,6,5,4),lwd=rep(2,4),legend=c(plot.legend[c(1,3,4,5)]),bty='n')
	}

