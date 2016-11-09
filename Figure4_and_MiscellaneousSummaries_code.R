
###set appropriate working directory for R###

setwd() 

###Functions for calculating various metrics for "success" of different approaches###

integ.bias <- function(dat.use,sigma.true){
	rect.width <- diff(seq(-15,10,length.out=500))[1]
	rect.area <- rect.width*abs(dat.use[,4])/sigma.true
	return(sum(rect.area))
}

integ.cov <- function(dat.use){
	rect.width <- diff(seq(-15,10,length.out=500))[1]
	rect.area <- rect.width*abs(dat.use[,9])
	return(sum(rect.area))
}

####creating matrices for storing summaries of various approaches###

res.mat <- matrix(nrow=16, ncol=4) #matrix to store results in
colnames(res.mat) <- c('Max ESSS','Max Bias per SD','Integrated Bias per SD','Integrated Coverage') #column names
rownames(res.mat) <- rep(c('MEM.pi_n','MEM.pi_e','comm.pr','shm.tau50'),4) #row names for table

####read in simulation results and calcluate summaries of performance###
####the simulation results used in the manuscript can be found in the Simulation Results sub-folder####
####simulation code that can be edited to accomodate new scenarios can be found in MEM_Simulations.R###

for(scenloop in 1:4){

	scen.plot <- scenloop
	prior.plot <- c('pi_n.txt','pi_e.txt')
	MEM.list <- list()
	for(i in 1:2){
		asd.MEM <- read.table(paste0('Simulation Results/MEM.estimators_scen',scen.plot,'_',prior.plot[i]),header=F)
		MEM.list[[i]] <- asd.MEM[order(asd.MEM[,1]),]
	}

	commpr <- read.table(paste0('Simulation Results/commprior.estimators_scen',scen.plot,'.txt'),header=F); commpr <- commpr[order(commpr[,1]),]
	meta.app <- read.table(paste0('Simulation Results/shm.estimators_scen',scen.plot,'.txt'),header=F); meta.app[which(meta.app[,6]<0),6] <- 0; meta.app <- meta.app[order(meta.app[,1]),]

	dat.list <- list(MEM.list[[1]],MEM.list[[2]],commpr,meta.app) #data to put in table

	for(i in 1:length(dat.list)){
		asd <- dat.list[[i]]
		asd <- asd[which(asd[,1]<=6),]
		res.mat[(scen.plot*length(dat.list) - length(dat.list)) + i,'Max ESSS'] <- max(asd[,6])
		res.mat[(scen.plot*length(dat.list) - length(dat.list)) + i,'Max Bias per SD'] <- max(abs(asd[,4])/4)
		res.mat[(scen.plot*length(dat.list) - length(dat.list)) + i,'Integrated Bias per SD'] <- integ.bias(dat.use=asd, sigma.true=4)
		res.mat[(scen.plot*length(dat.list) - length(dat.list)) + i,'Integrated Coverage'] <- integ.cov(dat.use=asd)
	}
}

####Create Figure 4###

####Sub-Figure (a) - max median ESSS as a function as a function of integrated bias per SD###

scen.col = c('black','gray65','red','blue')
plot.names <- c('MEM.pi_n','MEM.pi_e','comm.pr','shm.tau50') 

plot.legend <- c(expression(paste('MEM ',pi[n])),expression(paste('MEM ',pi[e])),'CP','SHM') #expression(paste('MEM ',pi[n^"l"])),  #expression(paste(plain(sin)*eta^2, posterior))
pch.list <- c(15,17,5,6) 

plot(x=-10,y=-10,ylim=c(0,175),xlim=c(0.05,.5),ylab='Max Median ESSS',xlab='Integrated Bias per SD')

x.col <- 3 #column in results to find variable on x-axis
y.col <- 1

for(i in 1:length(plot.names)){
	asd <- res.mat[which(rownames(res.mat)==plot.names[i]),]
	for(j in 1:4){
		points(x=asd[j,x.col],y=asd[j,y.col],col=scen.col[j],pch=pch.list[i],cex=1.5)
	}
}

legend('topright',pch=c(pch.list,15,15,15,15),col=c(rep('gray65',length(pch.list)),scen.col),legend=c(plot.legend,'Scenario 1','Scenario 2','Scenario 3','Scenario 4'),bty='y')
title('(a)')

####Sub-Figure (b) - max median ESSS as a function as a function of integrated coverage###

scen.col = c('black','gray65','red','blue')
plot.names <- c('MEM.pi_n','MEM.pi_e','comm.pr','shm.tau50') 
pch.list <- c(15,17,5,6) #16,...,0,1,2,

plot(x=-10,y=-10,ylim=c(0,175),xlim=c(0.93,0.95),ylab='Max Median ESSS',xlab='Integrated Coverage')

x.col <- 4 #column in res.mat to find variable on x-axis
y.col <- 1

for(i in 1:length(plot.names)){
	asd <- res.mat[which(rownames(res.mat)==plot.names[i]),]
	for(j in 1:4){
		points(x=asd[j,x.col]/(21),y=asd[j,y.col],col=scen.col[j],pch=pch.list[i],cex=1.5)
	}
}
title('(b)')

############################################
### Extra statistics for paper from res.mat referenced throughout the manuscript###

#calculate % of bias/ESSS for MEM vs. CP
mem.n <- res.mat[which(rownames(res.mat)=='MEM.pi_n'),]
mem.e <- res.mat[which(rownames(res.mat)=='MEM.pi_e'),]
cp <- res.mat[which(rownames(res.mat)=='comm.pr'),]

#maximum ESSS from all 4 scenarios
mem.n_max_esss <- max(mem.n[,'Max ESSS'])
mem.e_max_esss <- max(mem.e[,'Max ESSS'])
cp_max_esss <- max(cp[,'Max ESSS'])

#Percent ESSS for MEM compared to CP
mem.n_esss <- 100*(mem.n[,'Max ESSS'] - cp[,'Max ESSS']) / cp[,'Max ESSS']
mem.e_esss <- 100*(mem.e[,'Max ESSS'] - cp[,'Max ESSS']) / cp[,'Max ESSS']

#%bias decrease comparison between MEM and CP
mem.n_bias <- ( 100* (cp[,3] - mem.n[,3]) / cp[,3] )
mem.e_bias <- ( 100* (cp[,3] - mem.e[,3]) / cp[,3] )









