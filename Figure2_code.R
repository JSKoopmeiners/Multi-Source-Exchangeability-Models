
###set appropriate working directory for R###

setwd() 

###import MEM Functions used to create figures###

source('MEM_Functions.R')

######### Calculate weights to check for convergence######

ncur <- c(1,10,100,seq(1000,1000000,by=1000)) #these values result in spacing of model labels on figure to not overlap
svec = rep(4,4)
nvec = rep(100,3)

pr.use <- 'pi_n'

#####sub-figure (a): n for current study goes to infinity but n for supplementary sources are constant###
#####plot MEM model weights where model 2 appropriate with only primary sample n to inf to demonstrate need for all sample sizes to convergence to inf####

weight.tab.MEM(xvec=c(2,2,15,15),svec=svec,nvec=nvec,ncur=ncur,prior=pr.use,plot=T,cur.n_inf=T)
legend('topright',col=c('black',rainbow(7)),legend=c('Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7','Model 8'), lty=c(1,2,2,2,6,4,4,5),bty='n')
title('(a)')

#####sub-figure (b): n for current study and supplementary data sources go to infinity###
#####plot MEM model weights where model 2 is correct with all n's to inf to demonstrate consistency of MEM under pi_n with all n to inf#####

weight.tab.MEM(xvec=c(2,2,15,15),svec=svec,nvec=nvec,ncur=ncur,prior=pr.use,plot=T,cur.n_inf=F)
title('(b)')

#####sub-figure (c): n for current study and supplementary data sources go to infinity using priors on model-inclusion using standard BMA###
#####plot BMA model weights where model 5 is correct with all n's to inf to demonstrate lack of consistent convergence under pi_n#####

weight.tab.BMA(xvec=c(2,2,2,15),svec=svec,nvec=nvec,ncur=ncur,prior=pr.use,plot=T,cur.n_inf=F)
title('(c)')

#####sub-figure (d): n for current study and supplementary data sources using priors on source-inclusion probabilites###
#####plot MEM model weights where model 5 is correct with al n's to inf to demonstrate consistency of model weights#####

weight.tab.MEM(xvec=c(2,2,2,15),svec=svec,nvec=nvec,ncur=ncur,prior=pr.use,plot=T,cur.n_inf=F)
title('(d)')


