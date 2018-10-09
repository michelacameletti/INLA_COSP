remove(list=ls())

library(spdep)
library(rgeos) 
library(sp)
library(maptools) 
library(splancs)
library(INLA) 
library(fields)
library(RColorBrewer)
library(lattice)
library(parallel)
n.cores = detectCores()

cluster = T
run_parallel = T

if(!cluster) {
  setwd("~/Dropbox/Lavori condivisi/COSP_INLA/Simulations")
  #  setwd("E:/Dropbox/COSP_INLA")
  #  setwd("F:/dropbox/Dropbox/COSP_INLA")

  }else{
  
  setwd("/home/michela/COSP_INLA/Simulation_casestudy/Simulations")
  #setwd("/home/michela/COSP_INLA_Torino/Simulation_casestudy/Simulations")
  
  #  d = dirname(INLA:::inla.call.builtin())n_s 
  #  inla.setOption(inla.call = paste(d,"/inla.static", sep=""))
  #  inla.setOption(fmesher.call = paste(d,"/fmesher.static", sep=""))
}  


############################################################
# Load functions / data / scenario settings
############################################################
source(paste0(getwd(),"/R_Code_LATEST/Extra_functions.R"))
source(paste0(getwd(),"/R_Code_LATEST/Data_simulation.R"))
source(paste0(getwd(),"/R_Code_LATEST/Simulation_singlestep.R"))
scenario.settings = read.csv(paste0(getwd(),"/R_Code_LATEST/Scenarios_bhicv.csv"),header=T,sep=";")

#SET HERE the scenarios you want to run!
run_scenarios = seq(2,3) 

############################################################
# Set the parameters which do not change across simulations or scenarios
############################################################
n_simulations = 500 #n. of simulations for each scenario

n_random_fromposterior = 100 #how many values to sample from the parameter posterior distr of the exposure model
n_random_fromPPD = 200 #how many values to sample from the exposure posterior distribution (for each area)

highres_simulation_reggrid = 100 #(do not modify this!) it's the resolution of the grid of points - covering all the region - where we simulate the exposure

# Exposure model parameters:
b0 = 10 #intercept 
b1 = 2 #covariate coeff
sigma2_e =  0.1 #meas err variance
sigma2_xi = 0.5 #spatial variance

# Poisson model parameters:
gamma0_Poisson = NULL #intercept 
sigma2iid_Poisson = 0.05 #variance iid effect

#################################
# Run the simulation! 
#################################
if(run_parallel){
  
  # # Check if the file with the performance results exists (if yes, delete it)
  name_file_performance = paste0("Performance_scenarios_",range(run_scenarios)[1],"_",range(run_scenarios)[2],".csv")
  # if(file.exists(paste0(getwd(),"/Simulations/Output/",name_file_performance)))
  #   file.remove(paste0(getwd(),"/Simulations/Output/",name_file_performance))
  
  # Empty all the folders
  #unlink(paste0(getwd(),"/Data/*"))
  #unlink(paste0(getwd(),"/Figures/*"))
  #unlink(paste0(getwd(),"/Output/*"))
  
  mclapply(X=run_scenarios, 
           FUN=Simulation_singlestep,
           scenario.settings = scenario.settings,
           name_file_performance = name_file_performance,
           sigma2_xi = sigma2_xi,
           b0 = b0, b1 = b1,
           sigma2_e = sigma2_e,
           sigma2iid_Poisson= sigma2iid_Poisson,
           gamma0_Poisson = gamma0_Poisson,
           mc.cores=n.cores,
           mc.silent=FALSE)
  
} else {
  
  for(i in run_scenarios){ 
    # # Check if the file with the performance results exists in the Output folder (if YES, delete it)
    name_file_performance = paste0("Performance_scenarios_",range(run_scenarios)[1],"_",range(run_scenarios)[2],".csv")
    # if(file.exists(paste0(getwd(),"/Simulations/Output/",name_file_performance)))
    #   file.remove(paste0(getwd(),"/Simulations/Output/",name_file_performance))
    
    # Empty all the folders
    unlink(paste0(getwd(),"/Data/*"))
    unlink(paste0(getwd(),"/Figures/*"))
    unlink(paste0(getwd(),"/Output/*"))
    
    Simulation_singlestep(x=i, 
                          scenario.settings = scenario.settings,
                          name_file_performance = name_file_performance,
                          sigma2_xi = sigma2_xi,
                          b0 = b0, b1 = b1,
                          sigma2_e = sigma2_e,
                          sigma2iid_Poisson= sigma2iid_Poisson,
                          gamma0_Poisson = gamma0_Poisson)
  }
}



