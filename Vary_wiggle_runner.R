rm(list=ls())
CPU = system("hostname",intern=T)
.libPaths(c(.libPaths(),"~/R/"))
Args = commandArgs(trailingOnly = T)
debug = length(Args)==0

source('CRN_BO/Discrete_Optimizers.R')

# methods
method_names = c("UNAMED",
                 "KG",
                 "UNAMED",
                 "CRNKG",
                 "PWKG")

################################################################################

if(length(Args)>0){
  # Run from terminal

  cat(CPU,"\n\n", commandArgs(), "\n\n")
  cat(getwd(),"\n\n")
  
  method = as.integer(Args[1])
  rho    = seq(0, 1, 0.1)[as.integer(Args[2])]
  BOseed = as.integer(Args[3])
  
  
  myID = paste(c("WIG", Args), collapse = "_")

  Ns0    = 3
  NOISE  = 50^2


}else{
  # run from Rstudio

  if (CPU=="huanan") setwd("~/Dropbox/PhD/CRN/")
  else if (grepl("Book", CPU)) setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/")
  else setwd("/Users/pearce/CRN/")

  cat("Running locally \n")
  method = 5
  BOseed = 1
  myID   = NULL #"somemoretastyshit"
  Ns0    = 3
  rho    = 0.1
  NOISE  = 50^2
}


################################################################################
# Make test data
set.seed(BOseed)
TruePars = c(5, 100^2, 5, rho*NOISE, 0, 0.0001 + (1-rho)*NOISE)
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

AA = ALGORITHMS[[method]]$new(TestFun, ran, X_domain, Y_f, BOseed, myID)
AA$optimize(Budget0=5, Budget=100, hypers=TruePars)


