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
  cat(getwd(),"\n\n")
  # system("cp ../CRNMay2018.R ./", wait=T)
}


source('CRNMay2018.R')


#############################################################################################
#############################################################################################


# methods
# 1: iid uniform
# 2: iid kg
# 3: ?
# 4: CRN KG
# 5: CRN Frazier
method_names = c("UNI MLE",
                 "KG",
                 "UNI HP",
                 "CRNKG", 
                 "PWKG_grad_MLE", 
                 "PWKG_grad_HP",
                 "PWKG_MLE",
                 "PWKG_HP",
                 "CRNKG_MLE N0=100")



if(length(Args)>0){
  
  reps = 400
  Methods = rep(c(4, 5, 6, 7), each=reps)
  BOseeds  = rep(1:reps, len=length(Methods))
  Ns0 = 5
  
  myID = as.numeric(Args[1]) 
  
  method = Methods[myID]
  BOseed = BOseeds[myID]
  
}else{
  
  cat("Running locally \n")
  method = 5
  BOseed = 399
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
                    Timing= Timing
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
dims      = 8

# TestFun   = Build_Ambulance_Testfun(BOseed, numtestseeds=10000, runlength=1)[[1]]
TestFun   = Build_Xie_ATO_cpp_Testfun(BOseed, 2000, 1)[[1]]
XRAN      = matrix(c(0,20), 2, dims)


Budget0   = 20
Budget1   = 25

cat("TestFun Ambulance, method ", method_names[method], ", seed ", BOseed, "\n")

################################################################################################
################################################################################################


N0 = function()length(GP1$yd)
NS = function()length(unique(GP1$xd[,ncol(GP1$xd)]))

# Otherwise initialize the GP models

t0 = proc.time()[3]
XX   = UniformDesign_X(N0=Budget0, ran=XRAN, Ns=Ns0, TestFun=NULL, rounding=T, double=0) 
YY   = TestFun(XX)
eval_time = proc.time()[3] - t0

if(method%in%c(1, 2))  XX[,ncol(XX)] = 1:nrow(XX)

GP1  = CRNLHood$new(XX, YY, XRAN)

if(method<6){
  GP1$Refresh(learnHpars=1)
}else{
  GP1$Refresh(learnHpars=3)
}

fit_time = proc.time()[3] - t0 - eval_time

RecX = list()
RecX[[N0()]] = round(GP1$RecX())

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


while(length(GP1$yd)<Budget1){
  
  cat("Selecting point: ", length(GP1$yd)+1, "\n")
  if(method == 1){
    # iid uniform
    Budget = length(GP1$yd) + 1
    XX     = UniformDesign_X(Budget, XRAN, rounding=T, double=0)
    
    YY     = TestFun(XX)
    
    GP1    = CRNLHood$new(XX, YY, XRAN)
    
    GP1$Refresh()
    
  }else{
    t0 = proc.time()[3]
    
    # Sequential methods
    Xr = Build_ref_X(GP1, T)
    max_seed = max(GP1$xd[,GP1$dims+1])
    
    if(method==2){
      newx = MCMC_CRNKG_grad(list(GP1), check_Seeds=max_seed, Xr=Xr,
                             N0=1000, Na=5, maxevals=20)
      
    }else if(method==4 | method==6){
      
      check_seeds = sample(max_seed)[1:5]
      newx = MCMC_CRNKG_grad(list(GP1), Xr=Xr, check_Seeds=check_seeds,
                             N0=1000, Na=5, maxevals=20)
      
    }else if(method==5 | method==7){
      newx = MCMC_PWKG_grad(list(GP1), Xr=Xr,
                            N0=1000, Na=5, maxevals=20, 
                            PN0=4000, PNa=20, Pmaxevals=40)
      
    }
    
    KG_time = proc.time()[3] - t0
    
    newx = round(newx)
    
    newy = TestFun(newx)
    
    eval_time = proc.time()[3] - t0 - KG_time
    
    GP1$xd     = rbind(GP1$xd, newx)
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
    fit_time = proc.time()[3] - t0 - eval_time - KG_time
  }
  
  
  
  GP1$Lhood_Prep()
  Lhoods[[N0()]]       = GP1$Lhood_standard(GP1$HP)
  Hpar_History[[N0()]] = GP1$HP
  Pmeans[[N0()]]       = GP1$ymean
  RecX[[N0()]]         = round(GP1$RecX( oldrecx=tail(RecX,1)[[1]] ))
  Timing[[N0()]]       = c(KG_time, eval_time, fit_time)
  
  cat(" Recomended x: ", RecX[[N0()]], ", ")
  
  Cost[nrow(Cost)+1,] = c(N0(), TestFun( c(RecX[[N0()]], 0) ) )
  
  
  if(debug){plot(Cost); lines(Cost)}
  
  cat( "\n", method_names[method], " ", as.numeric(Cost[nrow(Cost),]), "\n")
  
  Checkpoint()
  
  cat("\n\n")
}


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
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/ATO_runner.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    # system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Ambulances.cpp godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    
  }
  NrunFileArray("/storage/maths/phrnaj/PhD/CRN/ATO_runner.R", 400*5)
  
  # NrestoreFileArray("Ambulances_runner.R", 1078, 200*5)
  
  ################################################
  # Orac
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/ATO_runner.R orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Ambulances.cpp orac:~/CRN/", wait = T)
  NOracFile("~/CRN/ATO_runner.R", 500)
  
  
  
  ################################################
  # Tinis
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R tinis:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/ATO_runner.R tinis:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Ambulances.cpp tinis:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp tinis:~/CRN/", wait = T)
  NTinisFile("~/CRN/ATO_runner.R", 500)
  
  
}

