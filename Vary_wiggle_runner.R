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
  # system("cp ../CRNMay2018.R ./", wait=T)
}

cat(getwd(),"\n\n")

source('CRNMay2018.R')
source('TestFuns.R')

# Rprof("tmpfile")
#############################################################################################
#############################################################################################


# methods
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
  # if(T){
  
  method = as.integer(Args[1])
  rho    = seq(0,1,0.1)[as.integer(Args[2])]
  BOseed = as.integer(Args[3])
  
  myID = paste(Args, collapse = "_")
  
  # reps = 800
  # Methods  = rep(rep(c(2, 4, 5),11), each=reps)
  # BOseeds  = rep(1:reps, len=length(Methods))
  # RHOS     = rep(seq(0,1,len=11), each=reps*3)
  # 
  # myID     = as.numeric(Args[1])
  # 
  # method = Methods[myID]
  # BOseed = BOseeds[myID]
  # rho    = RHOS[myID]
  
  Ns0    = 3
  NOISE  = 50^2
  TruePars = c(5, 100^2, 5, rho*NOISE, 0, (1-rho)*NOISE)
  
}else{
  
  cat("Running locally \n")
  method = 4
  BOseed = 1
  myID   = 1199
  Ns0    = 3
  rho    = 0
  TruePars = c(5, 100^2, 5, rho*50^2, 0, (1-rho)*50^2)
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
      cat("Done!\n")
    }
    
    SaveState = function(){
      cat("saving file ...", paste("EachData1/", myID, sep=""))
      Output = list(Cost   = Cost,
                    method = method, 
                    method_name = method_names[method],
                    TPars  = TruePars, 
                    seed   = BOseed, 
                    x      = GP1$xd, 
                    y      = GP1$yd,
                    maxY_f = max(Y_f),
                    rho    = rho,
                    TF     = "ATO Xie lx=10", 
                    info   = "PW with/o grad, all with dbl init, MLE hpars with Optimizer_2",
                    Hpars  = Hpar_History,
                    CPU    = system("hostname",intern=T),
                    myID   = myID,
                    .Random.seed = .Random.seed,
                    RecX   = RecX,
                    Lhoods = Lhoods,
                    Pmeans = Pmeans
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


cat("Method ", method, ", wiggle ", rho, ", seed ", BOseed, "\n")

####################################################
# Make test function
set.seed(BOseed)

X_domain = 1:100
SEkernel = function(x1,x2){
  TruePars[2]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[1]^2)
}
X_f = X_domain
Y_f = mvrnorm(1, rep(0,length(X_f)), SEkernel(X_f, X_f))

SEkernelw = function(x1,x2){
  TruePars[4]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[3]^2)
}
Offs_f = mvrnorm(60, rep(0, length(X_f)), SEkernelw(X_f,X_f))

TestFun = function(xs){
  
  
  TestFun_i = function(xs){
    eps = 0
    if(xs[2]>0) eps=Offs_f[xs[2], xs[1]] + rnorm(1,0,sqrt(TruePars[6]))
    Y_f[xs[1]] + eps
  }
  
  xs = matrix(xs, ncol=2)
  apply(xs,1,TestFun_i)
}

XRAN    = matrix(range(X_domain), 2)



#######################################################
# Iitialize the GP models

Budget0 = 5

Budget1 = 50


XX_init   = UniformDesign_X(N0=Budget0, ran=XRAN, Ns=3, TestFun=NULL, rounding=T, double=0) 
YY_init   = TestFun(XX_init)

if(method%in%c(1, 2))  XX_init[,2] = 1:nrow(XX_init)

GP1  = CRNLHood$new(XX_init, YY_init, XRAN)
N0 = function()length(GP1$yd)
NS = function()length(unique(GP1$xd[,ncol(GP1$xd)]))

GP1$Refresh(learnHpars=7, Hpars=TruePars)


RecX = list()
RecX[[N0()]] = X_domain[ which.max(GP1$MU(cbind(X_domain, 0))) ]

PPP  =  TestFun(c(RecX[[N0()]],0))
Cost =  data.frame(N=N0(), P= PPP)

Hpar_History = list()
Hpar_History[[N0()]] = GP1$HP

GP1$Lhood_Prep()
Lhoods = list()
Lhoods[[N0()]] = GP1$Lhood_standard(GP1$HP)

Pmeans = list()
Pmeans[[N0()]] = GP1$ymean


# Try and load old model, and save loaded/initial model
Checkpoint(tryload=T)

################################################################################################
################################################################################################
# Start Sequential allocation

cat("\n\nReady to Rock and Roll!\n\n")


while(length(GP1$yd)<Budget1){
  
  
  if(method == 1){
    
    # iid uniform
    Budget = length(GP1$yd) + 1
    XX     = UniformDesign_X(Budget, XRAN, rounding=T, double=0)
    
    YY     = TestFun(XX)
    
    GP1    = CRNLHood$new(XX, YY, XRAN)
    
    GP1$Refresh()
    
  }else{
    # Sequential methods
    Xr = cbind(X_domain, 0) #Build_ref_X(GP1, T)
    
    KG = Make_CRNKG_grad(GP1, Xr)
    PKG = Make_PWKG_grad(GP1, Xr)
    
    if(method==2){
      newx = X_domain[ which.max(sapply(X_domain, PKG$KG)) ]
      newx = c(newx, max(GP1$xd[,2])+1 )
      
    } else if(method==4){
      checkseeds = 1:(max(GP1$xd[,2])+1)
      KGvals = sapply(checkseeds, function(s)sapply(X_domain, function(xi)KG$KG(c(xi, s))))
      
      newx = X_domain[ which.max(apply(KGvals, 1, max)) ]
      news = which.max(apply(KGvals, 2, max))
      
      newx = c(newx,news)
      # browser()
      
    }else if(method==5){
      
      KG_vals = sapply(X_domain, PKG$KG)
      
      PW_vals = sapply(X_domain, function(x1)sapply(X_domain, function(x2)PKG$PWKG(c(x1,x2))))
      
      if(max(KG_vals) > max(PW_vals)){
        newx = X_domain[ which.max( KG_vals )]
        newx = c(newx, max(GP1$xd[,2])+1)
      }else{
        newx1 = X_domain[ which.max(apply(PW_vals,1,max)) ]
        newx2 = X_domain[ which.max(apply(PW_vals,2,max)) ]
        newx = cbind(c(newx1, newx2), max(GP1$xd[,2])+1 )
      }
      
      
    }
    
    newx = round(newx)
    
    newy = TestFun(newx)
    
    GP1$xd   = rbind(GP1$xd, newx)
    GP1$yd_o   = c(GP1$yd_o, newy)
    
    GP1$Refresh(learnHpars=7, Hpars=TruePars)
    
  }
  
  
  
  GP1$Lhood_Prep()
  Lhoods[[N0()]]       = GP1$Lhood_standard(GP1$HP)
  Hpar_History[[N0()]] = GP1$HP
  Pmeans[[N0()]]       = GP1$ymean
  RecX[[N0()]]         = X_domain[ which.max(GP1$MU(cbind(X_domain, 0))) ]
  
  cat(" Recomended x: ", RecX[[N0()]], ", ")
  
  PPP = TestFun( c(RecX[[N0()]], 0) )
  Cost[nrow(Cost)+1,] = c(N0(), PPP )
  
  
  if(debug){
    par(mfrow = c(2,1))
    X = 0:100
    MUY = GP1$MU(cbind(X,0))
    plot(X_f, Y_f, col="grey", lwd=2, type="l", ylab="Y", ylim=c(-200,200))
    lines(X, MUY, lwd=2)
    points(GP1$xd[,1], GP1$yd, col=GP1$xd[,2]+1, pch=19)
    
    # points(which.max(MUY)+1, max(MUY), pch=19, cex=3)
    
    plot(Cost); lines(Cost)
    
  }
  
  cat( "\n", method_names[method], " ", as.numeric(Cost[nrow(Cost),]), "\n")
  
  Checkpoint()
  
  cat("\n\n")
  # stop()
}

# Rprof()
# print(summaryRprof("tmpfile"))

cat("Finished and saved", myID, "###########################################\n")
cat("Finished and saved", myID, "###########################################\n")
cat("Finished and saved", myID, "###########################################\n\n\n")






















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
  if (grepl("huanan", CPU)){
    source("~/Dropbox/PhD/R/Nhtop.R")
    
    system("scp ~/Dropbox/PhD/CRN/CRNMay2018.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp ~/Dropbox/PhD/CRN/CRN_runner_GP_DiscComSph.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    
  }
  if (grepl("Book", CPU)){
    source("/Volumes/Datastorage/Dropbox/PhD/R/Nhtop.R")
    
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/Vary_rho_runner.R godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp godzilla:/storage/maths/phrnaj/PhD/CRN/", wait = T)
    
  }
  NrunFileArray("/storage/maths/phrnaj/PhD/CRN/Vary_rho_runner.R", 26400)
  
  
  # Orac
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/CRNMay2018.R orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/WSC_1D_runner.R orac:~/CRN/", wait = T)
  system("scp /Volumes/Datastorage/Dropbox/PhD/CRN/rcpp_ATO_kernel.cpp orac:~/CRN/", wait = T)
  NOracFile("~/CRN/WSC_1D_runner.R", 600)
  
  
}
