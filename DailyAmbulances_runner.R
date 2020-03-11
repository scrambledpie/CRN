rm(list=ls())
CPU = system("hostname",intern=T)
Args = commandArgs(trailingOnly = T)
debug = length(Args)==0
.libPaths(c(.libPaths(),"~/R/"))

if(!length(Args)%in%c(0,2))stop("Give a jobID and Budget>20 as arguments")



if(!debug){
  
  # Running from terminal
  cat(CPU,"\n\n", commandArgs(), "\n\n")


  cat("current dir: ", getwd(),"\n")

  # wd0 = getwd()
  # wd0_split = strsplit(wd0, "/")[[1]]
  # dirs = length(wd0_split) - 1
  # wd1 = paste(wd0_split[1:dirs], collapse="/")
  # setwd(wd1)
  # cat("changed dir to: ", getwd(),"\n\n")
  
  reps = 400
  Methods = rep(c(2, 4, 5, 6, 7), each=reps)
  BOseeds  = rep(1:reps, len=length(Methods))
  
  set.seed(1)
  JOBS = sample(length(Methods))
  myID = JOBS[as.numeric(Args[1])] 
  
  # Optimization Run
  method = Methods[myID]
  BOseed = BOseeds[myID]
  Ns0    = 5
  Budget = as.integer(Args[2])
  filename = paste(c(myID, "AMB", method, BOseed), collapse = "_")
  
}else{
  
  # running within Rstudio
  cat("Running locally \n")
  
  if (CPU=="huanan")        setwd("~/Dropbox/PhD/CRN/")
  if (grepl("Book", CPU))   setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/git_CRN/")
  if (grepl("pearce", CPU)) setwd("/Users/pearce/CRN/")
  
  # Optimization Run
  method   = 10
  BOseed   = 1
  Ns0      = 5
  Budget   = 500
  filename = NULL
}

# Make the TestFunction
source('TestFuns/TestFuns.R')
# TestFun = Build_Ambulance_Testfun(BOseed, numtestseeds=60000, runlength=1, NumCalls=5)[[1]]
TestFun = Build_DailyAmbulance_Testfun(BOseed, numtestseeds=60000, runlength=1, Simtime=1800)[[1]]
ran = attr(TestFun, 'ran')


# pick the algortihm from the list
source('CRN_BO/Optimizers.R')
source('CRN_BO/SOSA.R')
ALGORITHMS = c(1,
               BO_KG,
               1,
               BO_CRNKG_CS,
               BO_PWKG_CS,
               BO_CRNKG_CSW,
               BO_PWKG_CSW,
               BO_CRNKG_CS_allseeds,
               BO_CRNKG_CSW_allseeds,
               SOSA)

AA = ALGORITHMS[[method]]$new(TestFun, ran, BOseed, myID=filename, rounding=F)



# exectute the optimizer
AA$optimize(Budget0 = 20, Budget = Budget, num_ref_x=NULL, Ns=1,
            N0  = 1000, Na  = 5, maxevals  = 50,
            PN0 = 4000, PNa = 10, Pmaxevals = 100)
cat("Finished and Saved ", filename)


# 
# # Make the TestFunction
# source('TestFuns/TestFuns.R')
# TestFuns= Build_DailyAmbulance_Testfun(BOseed, numtestseeds=10000, runlength=1, Simtime=1800)
# 
# 
# cat(TestFuns[[1]](c(runif(6), 1)))
# 
