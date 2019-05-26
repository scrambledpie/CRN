rm(list=ls())
CPU = system("hostname",intern=T)

Args = commandArgs(trailingOnly = T)
debug = length(Args)==0

if(!debug)cat(CPU,"\n\n", commandArgs(), "\n\n")

.libPaths(c(.libPaths(),"~/R/"))


if(debug){
  if (CPU=="huanan") setwd("~/Dropbox/PhD/CRN/")
  else if (grepl("Book", CPU)) setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/")
  else setwd("/Users/pearce/CRN/")
}

cat(getwd(),"\n\n")


source('CRNMay2018.R')
source('TestFuns.R')

#############################################################################################
#############################################################################################
method_names = c("UNI MLE",
                 "KG",
                 "UNI HP",
                 "CRNKG-CS", 
                 "PWKG-CS", 
                 "CRNKG-CS-W",
                 "PWKG-CS-W",
                 "UNNAMED",
                 "UNNAMED")



if(length(Args)>0){
  
  reps     = 400
  Methods  = rep(c(2, 4, 5, 6, 7), each=reps)
  BOseeds  = rep(1:reps, len=length(Methods))
  Ns0      = 5
  
  set.seed(1)
  JOBS = sample(length(Methods))
  
  myID = JOBS[as.numeric(Args[1])]
  
  method = Methods[myID]
  BOseed = BOseeds[myID]
  
}else{
  
  cat("Running locally \n")
  method = 4
  BOseed = 1
  myID   = 1199
  Ns0    = 5
  
}

Checkpoint = function(tryload=F){
  if(length(Args)>0){
    
    cat("\nCheckpoint, ")
    
    LoadState = function(){
      cat("loading file ", paste("EachData1/", myID, sep=""), "....")
      Input = readRDS(paste("EachData1/", myID, sep=""))
      Cost <<- Input$Cost
      Hpar_History <<- Input$Hpars
      GP1 <<- CRNLHood$new(Input$x, Input$y, XRAN)
      GP1$Refresh(Hpars = Input$Hpars[[length(Input$y)]], learnHpars=7)
      BOseed <<- Input$seed
      method <<- Input$method
      .Random.seed <<- Input$.Random.seed
      RecX   <<- Input$RecX; if(is.null(RecX))RecX = list()
      Lhoods <<- Input$Lhoods
      Pmeans <<- Input$Pmeans
      Timing <<- Input$Timing
      cat("Done!\n")
    }
    
    SaveState = function(){
      cat("saving file ...", paste("EachData1/", myID, sep=""))
      Output = list(Cost  = Cost,
                    method = method, 
                    method_name = method_names[method],
                    seed  = BOseed, 
                    x     = GP1$xd, 
                    y     = GP1$yd, 
                    TF    = "ATO Xie lx=10", 
                    info  = "PW with/o grad, all with dbl init, MLE hpars with Optimizer_2",
                    Hpars = Hpar_History,
                    CPU   = system("hostname",intern=T),
                    myID  = myID,
                    .Random.seed = .Random.seed,
                    RecX  = RecX,
                    Lhoods= Lhoods,
                    Pmeans= Pmeans,
                    Timing=Timing
      )
      saveRDS(Output,paste("EachData1/", myID, sep=""))
      cat(" done\n")
    }
    
    if(tryload){
      if(file.exists(paste("EachData1/", myID, sep=""))){
        tryCatch(LoadState(),error=function(e){cat("Failed to load"); 1})
      }
    }
    
    SaveState()
    
  }
}



################################################################################################
################################################################################################
# Define Test Functions and their respetive input ranges

set.seed(BOseed)

TestFun   = Build_Ambulance_Testfun(BOseed, numtestseeds=10000, runlength=1)[[1]]
XRAN      = attr(TestFun, 'ran')
dims      = ncol(XRAN)


Budget0   = 20
Budget1   = 500

cat("TestFun Ambulance, method ", method_names[method], ", seed ", BOseed, "\n")

################################################################################################
################################################################################################
# Initialize the GP models

t1   = proc.time()[3]
XX   = UniformDesign_X(N0=Budget0, ran=XRAN, Ns=Ns0, TestFun=NULL, rounding=F, double=0) 
YY   = TestFun(XX)
eval_time = proc.time()[3] - t1

# If we are not using CRN then overwrite seeds
if(method%in%c(1, 2))  XX[,ncol(XX)] = 1:nrow(XX)

# make the model and learn the hypers depending on wiggle/CS
GP1  = CRNLHood$new(XX, YY, XRAN)
if(method<6){
  GP1$Refresh(learnHpars=1)
}else{
  GP1$Refresh(learnHpars=3)
}
fit_time = proc.time()[3] - t1 - eval_time

N0 = function()length(GP1$yd)
NS = function()length(unique(GP1$xd[,ncol(GP1$xd)]))

#########################################
# Initialize all the logs
RecX = list()
RecX[[N0()]] = GP1$RecX()

Cost = data.frame(N=N0(), P=TestFun(c(RecX[[N0()]],0)))

Hpar_History = list()
Hpar_History[[N0()]] = GP1$HP

GP1$Lhood_Prep()
Lhoods = list()
Lhoods[[N0()]] = GP1$Lhood_standard(GP1$HP)

Pmeans = list()
Pmeans[[N0()]] = GP1$ymean

Timing = list()
Timing[[N0()]] = c(0, eval_time, fit_time)

# Try and load old model, and save loaded/initial model
Checkpoint(tryload=T)

################################################################################################
################################################################################################
# Start Sequential allocation

cat("\n\nReady to Rock and Roll!\n\n")

# Rprof("tmpfile")

while(length(GP1$yd)<Budget1){
  
  
  if(method == 1){
    # iid uniform
    Budget  =  length(GP1$yd) + 1
    XX      =  UniformDesign_X(Budget, XRAN, rounding=F, double=0)
    
    YY      =  TestFun(XX)
    
    GP1     =  CRNLHood$new(XX, YY, XRAN)
    
    GP1$Refresh()
    
  }else{
    # Sequential methods
    
    t0 = proc.time()[3]
    Xr = Build_ref_X(GP1, T)
    new_seed = max(GP1$xd[,GP1$dims+1]) +1
    
    if(method==2){
      newx = MCMC_CRNKG_grad(list(GP1), check_Seeds=new_seed, Xr=Xr,
                             N0=2000, Na=10, maxevals=100)
      
    }else if(method==4 | method==6){
      check_seeds = sort(sample(new_seed)[1:min(new_seed, 5)])
      newx = MCMC_CRNKG_grad(list(GP1), Xr=Xr, check_Seeds=check_seeds,
                             N0=2000, Na=10, maxevals=100)
      
    }else if(method==5 | method==7){
      newx = MCMC_PWKG_grad(list(GP1), Xr=Xr,
                            N0=2000, Na=10, maxevals=100, 
                            PN0=8000, PNa=20, Pmaxevals=200)
      
    }
    KG_time = proc.time()[3]-t0
    
    newx = newx
    
    newy = TestFun(newx)
    
    eval_time = proc.time()[3] - t0 - KG_time
    
    GP1$xd   = rbind(GP1$xd, newx)
    GP1$yd_o   = c(GP1$yd_o, newy)
    
    OptimSteps = c(20:200, seq(205, 300, 5), seq(310, 400, 10), seq(420, 500, 20))
    
    if(N0()%in%OptimSteps){
      # Do a full hpar update
      if(method<6){
        GP1$Refresh(learnHpars=1)
      }else{
        GP1$Refresh(learnHpars=3)
      }
      
    }else{
      # Just finetune
      if(method<6){
        GP1$Refresh(learnHpars=2)
      }else{
        GP1$Refresh(learnHpars=4)
      }
    }
    
    fit_time = proc.time()[3] - t0 - KG_time - eval_time
    
  }
  
  
  # Update logs
  GP1$Lhood_Prep()
  Lhoods[[N0()]]       = GP1$Lhood_standard(GP1$HP)
  Hpar_History[[N0()]] = GP1$HP
  Pmeans[[N0()]]       = GP1$ymean
  RecX[[N0()]]         = GP1$RecX( oldrecx=tail(RecX,1)[[1]] )
  Timing[[N0()]]       = c(KG_time, eval_time, fit_time)
  
  # measure performance
  cat(" Recomended x: ", signif(RecX[[N0()]], 1), ", ")
  Cost[nrow(Cost)+1,] = c(N0(), TestFun( c(RecX[[N0()]], 0) ) )
  cat( "\n", method_names[method], " ", as.numeric(Cost[nrow(Cost),]), "\n")
  Checkpoint()
  
  
  if(debug){plot(Cost); lines(Cost)}
  
  
  
  cat("\n\n")
}

# Rprof()
# summaryRprof("tmpfile")

cat("Finished and saved", myID, "\n")






















if(F){
  source("/Volumes/DataStorage/Dropbox/PhD/R/Nhtop.R")
  system("scp /Volumes/DataStorage/Dropbox/PhD/CRN/CRNMay2018.R orac:~/CRN/",wait = T)
  
  if (grepl("huanan", CPU) | grepl("Book", CPU) ){
    source("~/Dropbox/PhD/R/Nhtop.R")
    system("scp ~/Dropbox/PhD/CRN/CRNMay2018.R orac:~/CRN/", wait = T)
  }
  NOracFile("~/CRN/CRNMay2018.R",1600)
  
  
  if (grepl("huanan", CPU)){
    source("~/Dropbox/PhD/R/Nhtop.R")
    
    system("scp ~/Dropbox/PhD/CRN/CRNMay2018.R tinis:~/CRN/", wait = T)
    system("scp ~/Dropbox/PhD/CRN/CRN_runner_GP_DiscComSph.R tinis:~/CRN/", wait = T)
    
    system("scp ~/Dropbox/PhD/CRN/CRNMay2018.R orac:~/CRN/", wait = T)
    system("scp ~/Dropbox/PhD/CRN/CRN_runner_GP_DiscComSph.R orac:~/CRN/", wait = T)
  }
  
  if (grepl("Book", CPU)){
    source("/Volumes/Datastorage/Dropbox/PhD/R/Nhtop.R")
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R tinis:~/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRN_runner_ATO.R tinis:~/CRN/", wait = T)
  }
  NTinisFile("~/CRN/CRN_runner_GP_DiscComSph.R", 1200)
  NOracFile("~/CRN/CRN_runner_GP_DiscComSph.R", 1200)
  
  
  
  
  ############################################
  ## CSC desktops
  if (grepl("Book", CPU)){
    source("/Volumes/Datastorage/Dropbox/PhD/R/Nhtop.R")
    
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Ambulances/Ambulances_runner.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Ambulances.cpp godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    
  }
  NrunFileArray("/storage/maths/phrnaj/PhD/CRN/Ambulances_runner.R", 400*5)
  
  NrestoreFileArray("Ambulances_runner.R", 1095, 400*5)
  
  
  # Orac
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Ruuner_8DGP_FixedHP.R orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp orac:~/CRN/", wait = T)
  NOracFile("~/CRN/Ruuner_8DGP_FixedHP.R", 2400)
  
  
}

