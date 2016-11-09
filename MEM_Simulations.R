nsim=NA

library(snowfall)
sfInit(cpus = 12, type = 'SOCK', parallel = T)

sfExportAll()

sfLibrary(gtools)
sfLibrary(matrixStats)

sfSource("MEM_Functions.R")

sfSapply(1:500, function(x) commensurateprior_calc_sim(x=c(x,-4,-4,-4),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen1'))

sfSapply(1:500, function(x) commensurateprior_calc_sim(x=c(x,-10,-10,2),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen2'))

sfSapply(1:500, function(x) commensurateprior_calc_sim(x=c(x,-10,-4,2),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen3'))

sfSapply(1:500, function(x) commensurateprior_calc_sim(x=c(x,-10,-9.25,2),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen4'))

sfSapply(1:500, function(x) shm_calc_sim(x=c(x,-4,-4,-4),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen1'))

sfSapply(1:500, function(x) shm_calc_sim(x=c(x,-10,-10,2),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen2'))

sfSapply(1:500, function(x) shm_calc_sim(x=c(x,-10,-4,2),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen3'))

sfSapply(1:500, function(x) shm_calc_sim(x=c(x,-10,-9.25,2),S=c(4,4,4,4),N=c(100,100,100,100),scen='scen4'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-4,-4,-4),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_n',scen='scen1'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-10,-10,2),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_n',scen='scen2'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-10,-4,2),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_n',scen='scen3'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-10,-9.25,2),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_n',scen='scen4'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-4,-4,-4),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_e',scen='scen1'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-10,-10,2),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_e',scen='scen2'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-10,-4,2),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_e',scen='scen3'))

sfSapply(1:500, function(x) MEM_sim_calc(x=c(x,-10,-9.25,2),S=c(4,4,4,4),N=c(100,100,100,100),prior='pi_e',scen='scen4'))

sfStop()