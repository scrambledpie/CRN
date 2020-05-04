library(R6)

# setwd("/home/michael/CRN/")

# setwd("/mnt/huanan/home/michael/CRN/")
setwd("/home/michael/CRN/")


source("CRN_BO/GP_model.R")

# First load up the test functions
source("TestFuns/TestFuns.R")

# Make the TestFunction
TESTFUNS = c(Build_Xie_ATO_cpp_Testfun,
             Build_Ambulance_Testfun,
             Build_DailyAmbulance_Testfun,
             Build_DailyAmbulanceSum_Testfun)

BOseed = 1
problem = 4

TestFun = TESTFUNS[[problem]](BOseed, numtestseeds=2000, runlength=1)[[1]]
ran = attr(TestFun, 'ran') # bounds of the search space.

Cross_validate = function(GP){
  # keeps the hyperparameters, removes one data points at a time and 
  # compute the delta to the missing point.
  
  n = length(GP$yd)
  
  full_X = GP$xd
  full_Y = GP$yd
  
  # browser()
  cv_delta_Z = sapply(1:n, function(i){
    # modify data and refit GP
    GP$xd = full_X[-i,]
    GP$yd_o = full_Y[-i]
    
    # GP$xd
    GP$Refresh(learnHpars=0)
    
    # compute z-score of missing point
    cv_X = full_X[i, , drop=F]
    cv_Y = full_Y[i]
    
    MU = GP$MU(cv_X)[1, 1]
    SIG = GP$COV(cv_X, cv_X)[1, 1]
    
    cv_delta = cv_Y - MU
    cv_Z = cv_delta / sqrt(SIG)
    
    return(c(cv_delta, cv_Z))
  })
  
  return(cv_delta_Z)
}


# First make a dataset of points
n_points = 200


for(seed in 1:5){
  set.seed(seed)
  
  X0 = LHCran_seeds(N0=n_points, ran=ran, ncats=5)
  X0 = round(X0)
  Y0 = TestFun(X0)
  cat("\nGot the dataset\n\n")
  
  
  
  # learnHpars dictates what type of GP model we learn
  # 0: keep previous hyperparameters
  # 1: Compound Shperic, full optimisation
  # 2: compound Spheric, just fine tune update hpars
  # 3: wiggles
  # 4: wiggles fine tune
  # 5: unique seeds i.i.d. noise GP full optim
  # 6: unique seeds i.i.d. noise GP finetune
  # 7: use given hyperparameter, "Hpars"
  
  
  # IID MODEL FIRST
  IA_GP = CRNLHood$new(X0, Y0, ran)
  IA_GP$Refresh(learnHpars = 5); cat("\nMade the IID model\n\n")
  IA_CV = Cross_validate(IA_GP); cat("got the CV for IID\n\n\n")
  
  
  # OFFSETS MODEL 
  OA_GP = CRNLHood$new(X0, Y0, ran)
  OA_GP$Refresh(learnHpars = 1); cat("\nMade the OFF model\n\n")
  OA_CV = Cross_validate(OA_GP); cat("got the CV for OFF\n\n\n")
  
  
  # OFFSETS + WIGGLES MODEL
  WA_GP = CRNLHood$new(X0, Y0, ran)
  WA_GP$Refresh(learnHpars = 3); cat("\nMade the OFF+WIG model\n\n")
  WA_CV = Cross_validate(WA_GP); cat("got the CV for OFF WIG\n\n\n")
  
  output = list(IA_CV=IA_CV, OA_CV=OA_CV, WA_CV=WA_CV)
  saveRDS(output, paste("paper_plots/CV_AMB_SUM/CV_", seed, sep=""))
  cat("\n##########################################\n###############  SAVED SEED ", seed, "###########\n\n\n\n")
  
  nr = min(ncol(IA_CV), ncol(OA_CV), ncol(WA_CV))
  
  delta_vals = cbind(IA_CV[1, 1:nr], OA_CV[1, 1:nr], WA_CV[1, 1:nr])
  z_vals = cbind(IA_CV[2, 1:nr], OA_CV[2, 1:nr], WA_CV[2, 1:nr])
  
  par(mfrow=c(1, 4))
  # boxplot(abs(delta_vals), main="AMBULANCE, 200 LHC points, 5 seeds, abs( delta )")
  boxplot(abs(delta_vals), main="ATO, 200 LHC points, 5 seeds, abs( delta )")
  
  measures = cbind(
    apply(abs(delta_vals), 2, mean),
    apply(abs(delta_vals), 2, median),
    apply(abs(delta_vals)^2, 2, mean),
    apply(abs(delta_vals)^2, 2, median)
  )
  
  # model_names = c("AMB, IID", "AMB, OFF", "AMB, OFF+WIG")
  model_names = c("ATO, IID", "ATO, OFF", "ATO, OFF+WIG")
  print(measures)
  for (i in 1:3){
    qqnorm(z_vals[, i], asp=1, main=model_names[i])
    lines(c(-10, 10), c(-10, 10))
  }
  
  
  

}


