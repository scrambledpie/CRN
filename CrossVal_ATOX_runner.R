# This code opens up each of the old files and runs a cross validation experiment on the
# final data and parameters.



rm(list=ls())
CPU = system("hostname",intern=T)
Args = commandArgs(trailingOnly = T)
debug = length(Args)==0
.libPaths(c(.libPaths(),"~/R/"))

if(!length(Args)%in%c(0,2))stop("Give a jobID and Budget>20 as arguments")


library(R6)
source("CRN_BO/GP_model.R")

# Running from terminal
cat(CPU,"\n\n", commandArgs(), "\n\n")
cat("current dir: ", getwd(),"\n")

# parse arguments
dirname = Args[1]
myID = as.numeric(Args[2]) + 1


filename = paste(c(myID, "CV_ATO"), collapse = "_")




Cross_validate = function(GP){
  # keeps the hyperparameters, removes one data points at a time and 
  # compute the delta to the missing point.
  
  n = 5 #length(GP$yd)
  
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

restore_data_perform_CV = function(data){
    # reload data and hypers from the data
    # m sets the noise model.

    X0 = data$x
    Y0 = data$y
    HP0 = data$Hpars[[length(data$Hpars)]]
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
    HP_IID = HP0
    HP_IID[18] = 0 #  wiggles
    HP_IID[19] = 0 #  offset
    HP_IID[20] = HP0[18] + HP0[19] + HP0[20] # white noise
    IA_GP = CRNLHood$new(X0, Y0, ran)
    IA_GP$Refresh(learnHpars = 7, Hpars=HP_IID); cat("\nMade the IID model\n\n")
    IA_CV = Cross_validate(IA_GP); cat("got the CV for IID\n\n\n")


    # OFFSETS MODEL
    HP_OFF = HP0
    HP_OFF[18] = 0 #  wiggles
    HP_OFF[19] = 0.5*HP0[18] + HP0[19] #  offset
    HP_OFF[20] = 0.5*HP0[18] + HP0[20] # white noise
    OA_GP = CRNLHood$new(X0, Y0, ran)
    OA_GP$HP = HP_OFF
    OA_GP$Refresh(learnHpars = 1); cat("\nMade the OFF model\n\n")
    OA_CV = Cross_validate(OA_GP); cat("got the CV for OFF\n\n\n")


    # OFFSETS + WIGGLES MODEL
    WA_GP = CRNLHood$new(X0, Y0, ran)
    WA_GP$Refresh(learnHpars = 7, Hpars=HP0); cat("\nMade the OFF+WIG model\n\n")
    WA_CV = Cross_validate(WA_GP); cat("got the CV for OFF WIG\n\n\n")

    return(list(IID_CV=IA_CV, OFF_CV=OA_CV, WIG_CV=WA_CV))

}



# restore results data from this run
res_dir = "/home/maths/phrnaj/RESULTS/CRN/181.phrnaj.2020-04-27.22:36:24.265b316/res/"
res_files = dir(res_dir)
res_data = readRDS(res_files[myID])

CV_results = restore_data_perform_CV(res_Data)

saveRDS(CV_results, filename)
cat("\n Finished and Saved", myID)

