rm(list=ls())
CPU = system("hostname",intern=T)
.libPaths(c(.libPaths(),"~/R/"))
Args = commandArgs(trailingOnly = T)
debug = length(Args)==0



# methods
method_names = c("UNAMED",
                 "KG",
                 "UNAMED",
                 "CRNKG",
                 "PWKG")

################################################################################

if(length(Args)==3){
  # Run from terminal
  
  cat(CPU,"\n\n", commandArgs(), "\n\n")
  cat(getwd(),"\n\n")
  cat("\nRunning with 3 arguments\n")
  
  method = as.integer(Args[1])
  rho    = seq(0, 1, 0.1)[as.integer(Args[2])]
  BOseed = as.integer(Args[3])
  
  
  myID = paste(c("OFF", Args), collapse = "_")
  
  cat(paste("myID:", myID))
  
  Ns0    = 3
  NOISE  = 50^2

} else if(length(Args)==1){
  # Run from terminal on CSC machines
  
  cat(CPU,"\n\n", commandArgs(), "\n\n")
  cat(getwd(),"\n\n")
  cat("\nRunning with 1 argument\n")
  
  
  reps = 100
  rhosteps = 11
  
  METHODS = rep(c(2,4,5), each=reps*rhosteps)
  RHO     = rep(rep(seq(0,1,len=rhosteps), each=reps), 3)
  SEEDS   = rep(rep(1:reps, rhosteps), 3)
  
  XP_id = as.integer(Args[1])
  
  method = METHODS[XP_id]
  rho    = RHO[XP_id]
  BOseed = SEEDS[XP_id]
  
  
  myID = paste(c("OFF", method, rho, BOseed), collapse = "_")
  
  cat(paste("myID:", myID))
  
  Ns0    = 3
  NOISE  = 50^2

  if( !any("CRN_BO"==dir()) ){
    cat("\ncopying CRN_BO"); system("cp -r ../CRN_BO ./") 
  }
  Sys.sleep(0.25)
  

}else{
  # run from Rstudio
  cat("\nRunning with no arguments\n")
  
  if (CPU=="huanan") setwd("~/Dropbox/PhD/CRN/CRN/")
  else if (grepl("Book", CPU)) setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/")
  else setwd("/Users/pearce/CRN/")
  
  cat("Running locally \n")
  method = 4
  BOseed = 97
  myID   = NULL #"somemoretastyshit"
  Ns0    = 3
  rho    = 0
  NOISE  = 50^2
}



source('CRN_BO/Discrete_Optimizers.R')

# run 10 experiments for each run
for (i in BOseed + 1000*(0:9)){
  
  myID_i = paste(myID, i, sep="_")
  
  ################################################################################
  # Make test data
  set.seed(i)
  TruePars = c(10, 100^2, 5, 0, rho*NOISE, 0.0001 + (1-rho)*NOISE)
  Nseeds   = 115
  X_domain = 1:100
  
  SEkernel = function(x1,x2) TruePars[2]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[1]^2)
  X_f = X_domain
  Y_f = mvrnorm(1, rep(0,length(X_f)), SEkernel(X_f, X_f))
  Y_f = matrix(Y_f, 100, Nseeds+1)
  
  SEkernelw = function(x1,x2)TruePars[4]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[3]^2)
  Y_wiggles = t(mvrnorm(Nseeds, rep(0, length(X_f)), SEkernelw(X_f,X_f)))
  Y_wiggles = cbind(0, Y_wiggles)
  
  Y_Offsets = c(0, rnorm(Nseeds, 0, sqrt(TruePars[5])))
  Y_Offsets = matrix(Y_Offsets , 100, Nseeds+1, byrow = T)
  
  Y_white   = matrix(rnorm(100*Nseeds, 0, sqrt(TruePars[6])), 100, Nseeds)
  Y_white   = cbind(0, Y_white)
  
  TestFunY  = Y_f + Y_wiggles + Y_Offsets + Y_white
  
  TestFun = function(xs){
    if (any(xs != round(xs))) stop("give integer inputs")
    xs = matrix(xs, ncol=2)
    apply(xs,1,function(xs) TestFunY[xs[1], xs[2]+1])
  }
  ran    = matrix(range(X_domain), 2)
  
  
  
  ################################################################################
  # Run the optimizers
  
  ALGORITHMS = c(1,
                 BO_KG_DISC,
                 1,
                 BO_CRNKG_DISC,
                 BO_PWKG_DISC)
  # browser()
  AA = ALGORITHMS[[method]]$new(TestFun, ran, X_domain, Y_f, i, myID_i)
  AA$optimize(Budget0=5, Budget=50, hypers=TruePars)
  
}












if(F){
  
  
  if (CPU=="huanan"){
    system("rsync -r -v /home/michael/Dropbox/PhD/CRN/CRN/ godzilla:/storage/maths/phrnaj/PhD/CRN/")
  }
  
  NrunFileArray("/storage/maths/phrnaj/PhD/CRN/Vary_rho_runner.R", 3300)
  
  
  
}