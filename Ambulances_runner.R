rm(list=ls())
CPU = system("hostname",intern=T)
Args = commandArgs(trailingOnly = T)
debug = length(Args)==0
.libPaths(c(.libPaths(),"~/R/"))

if(!length(Args)%in%c(0,2))stop("Give a jobID and Budget>20 as arguments")

if(!debug){
  # Running from terminal
  cat(CPU,"\n\n", commandArgs(), "\n\n")
  cat(getwd(),"\n\n")
  
  reps = 400
  Methods = rep(c(2, 8, 9, 6, 7), each=reps)
  BOseeds  = rep(1:reps, len=length(Methods))
  
  set.seed(1)
  JOBS = sample(length(Methods))
  myID = JOBS[as.numeric(Args[1])] 
  
  # Optimization Run
  method = Methods[myID]
  BOseed = BOseeds[myID]
  Ns0 = 5
  Budget = as.integer(Args[2])
  filename = paste(c(myID, "AMB", method, BOseed), collapse = "_")
  
}else{
  # running within Rstudio
  cat("Running locally \n")
  
  if (CPU=="huanan")        setwd("~/Dropbox/PhD/CRN/")
  if (grepl("Book", CPU))   setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/git_CRN/")
  if (grepl("pearce", CPU)) setwd("/Users/pearce/CRN/")
  
  # Optimization Run
  method   = 2
  BOseed   = 1
  Ns0      = 5
  Budget   = 25
  filename = NULL
}

# Make the TestFunction
source('TestFuns/TestFuns.R')
TestFun = Build_Ambulance_Testfun(1, 10000, 1)[[1]]
ran = attr(TestFun, 'ran')


# pick the algortihm from the list
source('CRN_BO/Optimizers.R')
ALGORITHMS = c(1, 
               BO_KG, 
               1, 
               BO_CRNKG_CS, 
               BO_PWKG_CS, 
               BO_CRNKG_CSW,
               BO_PWKG_CSW,
               BO_CRNKG_CS_allseeds,
               BO_CRNKG_CSW_allseeds)


# exectute the optimizer
AA = ALGORITHMS[[method]]$new(TestFun, ran, BOseed, myID=filename, rounding=F)
AA$optimize(Budget0 = 20, Budget = Budget,
            N0  = 1000, Na  = 10, maxevals  = 100,
            PN0 = 4000, PNa = 20, Pmaxevals = 200)
cat("Finished and Saved ", filename)


