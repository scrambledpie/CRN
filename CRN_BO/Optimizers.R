library(R6)
source("CRN_BO/GP_model.R")
source("CRN_BO/Acquisition_funs.R")
######################################################################
# A generic parent class for any optimizer to be used
OptimizerBase = R6Class("OptBase",
  public=list(
    TestFun    = NULL,
    ran        = NULL,
    X_domain   = NULL,
    RecX       = list(),
    Cost       = data.frame(N=numeric(0), P=numeric(0)),
    rounding   = FALSE,
    initialize = function(TestFun, ran, rounding=F, X_domain=NULL){
      self$TestFun = TestFun
      self$ran = ran
      self$rounding = rounding
      self$X_domain = X_domain
    }
  )
)


######################################################################
# A generic BO base class to be inherited by all BO methods
BO_base = R6Class("BO_base",inherit = OptimizerBase,
  public=list(
    Hpars  = list(),
    Lhoods = list(), 
    Ymean  = list(),
    Timing = list(),
    GP     = NULL,
    myID   = NULL,
    BOseed = NULL,
    method = NULL,
    Ns0    = NULL,
    initialize = function(testfun, ran, BOseed=NULL, 
                          myID=NULL, Ns0=5,
                          rounding=F
                          ){
      super$initialize(testfun, ran, rounding)
      self$Ns0     = Ns0
      self$BOseed  = BOseed
      self$myID    = myID
      set.seed(BOseed)
    },
    UpdateLogs = function(KG_time=0, eval_time=0, fit_time=0){
      N = length(self$GP$yd)
      
      self$GP$Lhood_Prep()
      self$Lhoods[[N]] = self$GP$Lhood_standard(self$GP$HP)
      
      self$Hpars[[N]]  = self$GP$HP
      self$Ymean[[N]]  = self$GP$ymean
      self$Timing[[N]] = c(KG_time, eval_time, fit_time)
      
      if(length(self$RecX)>0){
        oldxr = tail(self$RecX, 1)[[1]]
        self$RecX[[N]]   = self$GP$RecX(oldrecx=oldxr)
      }else{
        self$RecX[[N]]   = self$GP$RecX()
      }
      
      Rx               = matrix(c(self$RecX[[N]], 0), 1)
      NCi              = nrow(self$Cost)+1
      self$Cost[NCi,]  = c(N, self$TestFun(Rx))
    },
    Checkpoint = function(tryload=F){
      loadsuccess = F
      if(!is.null(self$myID)){
        
        myID = self$myID
        
        cat("\nCheckpoint, ")
        
        LoadState = function(){
          # get the file
          cat("loading file ", paste("EachData1/", myID, sep=""), "....")
          Input = readRDS(paste("EachData1/", myID, sep=""))
          
          # now copy everything over!
          self$Cost    <<- Input$Cost
          self$GP      <<- CRNLHood$new(Input$x, Input$y, self$ran)
          self$GP$Refresh(learnHpars=7, Hpars = Input$Hpars[[length(Input$y)]])
          self$BOseed  <<- Input$seed
          self$method  <<- Input$method
          .Random.seed <<- Input$.Random.seed
          self$RecX    <<- Input$RecX
          self$Lhoods  <<- Input$Lhoods
          self$Ymean   <<- Input$Ymean
          self$Timing  <<- Input$Timing
          self$Hpars   <<- Input$Hpars
          
          cat("Done!\n")
          
          return(T)
        }
        
        SaveState = function(){
          
          cat("saving file ...", paste("EachData1/", myID, sep=""))
          Output = list(CPU    = system("hostname",intern=T),
                        Cost   = self$Cost,
                        method = self$method, 
                        seed   = self$BOseed, 
                        x      = self$GP$xd, 
                        y      = self$GP$yd,
                        Hpars  = self$Hpars,
                        myID   = self$myID,
                        RecX   = self$RecX,
                        Lhoods = self$Lhoods,
                        Ymean  = self$Ymean,
                        Timing = self$Timing,
                        .Random.seed = .Random.seed,
                        methodname = class(self)[1]
          )
          saveRDS(Output,paste("EachData1/", myID, sep=""))
          cat(" done\n")
        }
        
        
        if(tryload){
          if(file.exists(paste("EachData1/", myID, sep=""))){
            loadsuccess = T
            tryCatch(LoadState(),error=function(e){cat("Failed to load"); loadsuccess=F; return(F)})
          }
        }else{
          SaveState()
        }
      }
      return(loadsuccess)
    },
    base_optimize = function(Budget0=20, Budget=500, get_next_x=NULL, learn_kernel=NULL, Ns0=5){
      
      
      if(is.null(get_next_x))stop("provide a get_next_x() function!")
      if(is.null(learn_kernel))stop("which noise kernel do we use? 0:iid, 1:CS,  2:CS+wig")
      
      # for each kernel a different index is passed to the hyper learning
      if(learn_kernel==0) kk = 5
      if(learn_kernel==1) kk = 1
      if(learn_kernel==2) kk = 3
      OptimSteps = c(20:200, seq(205, 300, 5), seq(310, 400, 10), seq(420, 500, 20))
      
      loadsuccess = self$Checkpoint(tryload = T)
      
      
      # if we failed to load, then initialize
      if(!loadsuccess){
        cat("not loaded. ")
        t0       = proc.time()[3]
        
        cat("\nSelecting points 1,...,", Budget0, "\n")
        # get first initilialization points
        X_init   = UniformDesign_X(N0=Budget0, ran=self$ran, Ns=Ns0, 
                                   TestFun=NULL, rounding=self$rounding, double=0) 
        KG_time = proc.time()[3] - t0
        
        # get objective function values
        Y_init    = self$TestFun(X_init)
        eval_time = proc.time()[3] - t0 - KG_time
        
        # fit the model
        self$GP = CRNLHood$new(X_init, Y_init, self$ran)
        self$GP$Refresh(learnHpars = kk)
        fit_time = proc.time()[3] - t0 - KG_time - eval_time
        
        
        self$UpdateLogs(KG_time, eval_time, fit_time)
        
        RX = signif(tail(self$RecX,1)[[1]], 3)
        cat("Logs updated, RecX:", RX,", performance:", as.numeric(tail(self$Cost,1)), "\n\n\n")
      }
        
      # Now do the sequential bit!
      while (length(self$GP$yd)<Budget){

        N =  length(self$GP$yd)+1
        cat("Selecting point: ",N, "\n")
        t0 = proc.time()[3]
        
        # Get the new x to evaluate
        newx = get_next_x()
        if(self$rounding) newx = round(newx)
        KG_time = proc.time()[3] - t0
        
        
        # Evaluate the objective function
        newy = self$TestFun(newx)
        eval_time = proc.time()[3] - t0 - KG_time
        
        
        # Update the data and either do a full hpar update or just finetune
        self$GP$xd      = rbind(self$GP$xd, newx)
        self$GP$yd_o   = c(self$GP$yd_o, newy)
        finetune        = !(N%in%OptimSteps)
        self$GP$Refresh(learnHpars = kk + finetune)
        fit_time        = proc.time()[3] - t0 - eval_time - KG_time
        
        self$UpdateLogs(KG_time, eval_time, fit_time)
        
        self$Checkpoint()
        
        RX = signif(tail(self$RecX,1)[[1]], 3)
        cat("Logs updated, RecX:", RX,", performance:", as.numeric(tail(self$Cost,1)), "\n\n\n")
      }
      
    }    
  )
)
######################################################################
######################################################################
######################################################################
######################################################################
# Each BO method

BO_KG = R6Class("BO_KG",
  inherit = BO_base,
  public = list(
    method = 2,
    optimize = function(Budget0=20, Budget=500, N0=1000, Na=10, maxevals=100, num_ref_x=NULL, ...){
     
      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        Xr = Build_ref_X(self$GP, T, num_ref_x)
        check_seeds = max(self$GP$xd[,self$GP$dims+1]) + 1
        newx = MCMC_CRNKG_grad(list(self$GP), Xr=Xr, check_Seeds=check_seeds,
                               N0=N0, Na=Na, maxevals=maxevals)
        return(newx)
      }
      
      # call the optimizer with the suggestion fuction and the kernel
      self$base_optimize(Budget0, Budget, get_next_x, 0)
    }
  )
)

######################################################################

BO_CRNKG_CS = R6Class("BO_CRNKG_CS",
  inherit = BO_base,
  public = list(
    method = 8,
    optimize = function(Budget0=20, Budget=500, N0=1000, Na=10, maxevals=100,num_ref_x=NULL,Ns=NULL,...){
      
      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        Xr = Build_ref_X(self$GP, T, num_ref_x)
        new_seed = max(self$GP$xd[,self$GP$dims+1]) + 1
        check_seeds = sample(new_seed)[1:min(Ns, new_seed)]
        newx = MCMC_CRNKG_grad_cheap(list(self$GP), Xr=Xr, check_Seeds=check_seeds,
                               N0=N0, Na=Na, maxevals=maxevals)
        return(newx)
      }
      
      # call the optimizer with the suggestion fuction and the kernel
      self$base_optimize(Budget0, Budget, get_next_x, 1)
    }
  )
)

######################################################################

BO_CRNKG_CSW = R6Class("BO_CRNKG_CSW",
  inherit = BO_base,
  public = list(
    method = 9,
    optimize = function(Budget0=20, Budget=500, N0=1000, Na=10, maxevals=100,num_ref_x=NULL, Ns=NULL,...){
      
      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        Xr = Build_ref_X(self$GP, T, num_ref_x)
        new_seed = max(self$GP$xd[,self$GP$dims+1]) + 1
        check_seeds = sample(new_seed)[1:min(Ns, new_seed)]
        newx = MCMC_CRNKG_grad_cheap(list(self$GP), Xr=Xr, check_Seeds=check_seeds,
                               N0=N0, Na=Na, maxevals=maxevals)
        return(newx)
      }
      
      # call the optimizer with the suggestion fuction and the kernel
      self$base_optimize(Budget0, Budget, get_next_x, 2)
    }
  )
)

######################################################################

BO_CRNKG_CS_allseeds = R6Class("BO_CRNKG_CS_allseeds",
  inherit = BO_base,
  public = list(
    method = 4,
    optimize = function(Budget0=20, Budget=500, N0=1000, Na=10, maxevals=100,num_ref_x=NULL, Ns=NULL,...){
      
      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        Xr = Build_ref_X(self$GP, T, num_ref_x)
        new_seed = max(self$GP$xd[,self$GP$dims+1]) + 1
        check_seeds = 1:new_seed
        newx = MCMC_CRNKG_grad(list(self$GP), Xr=Xr, check_Seeds=check_seeds,
                               N0=N0, Na=Na, maxevals=maxevals)
        return(newx)
      }
      
      # call the optimizer with the suggestion fuction and the kernel
      self$base_optimize(Budget0, Budget, get_next_x, 1)
    }
  )
)

######################################################################

BO_CRNKG_CSW_allseeds = R6Class("BO_CRNKG_CSW_allseeds",
  inherit = BO_base,
  public = list(
   method = 6,
   optimize = function(Budget0=20, Budget=500, N0=1000, Na=10, maxevals=100,num_ref_x=NULL,Ns=NULL,...){
     
     # define a function that suggests the next x to evaluate
     get_next_x = function(){
       Xr = Build_ref_X(self$GP, T, num_ref_x)
       new_seed = max(self$GP$xd[,self$GP$dims+1]) + 1
       check_seeds = 1:new_seed
       newx = MCMC_CRNKG_grad(list(self$GP), Xr=Xr, check_Seeds=check_seeds,
                              N0=N0, Na=Na, maxevals=maxevals)
       return(newx)
     }
     
     # call the optimizer with the suggestion fuction and the kernel
     self$base_optimize(Budget0, Budget, get_next_x, 2)
   }
  )
)
######################################################################

BO_PWKG_CS = R6Class("BO_PWKG_CS",
  inherit = BO_base,
  public = list(
    method = 5,
    optimize = function(Budget0=20, Budget=500,num_ref_x=NULL, 
                        N0=1000, Na=10, maxevals=100,
                        PN0=4000, PNa=40, Pmaxevals=200,...){
     
      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        Xr = Build_ref_X(self$GP, T, num_ref_x)
        newx = MCMC_PWKG_grad(list(self$GP), Xr=Xr,
                               N0=N0, Na=Na, maxevals=maxevals, 
                               PN0=PN0, PNa=PNa, Pmaxevals=Pmaxevals)
        return(newx)
      }
     
     # call the optimizer with the suggestion fuction and the kernel
     self$base_optimize(Budget0, Budget, get_next_x, 1)
   }
  )
)

######################################################################

BO_PWKG_CSW = R6Class("BO_PWKG_CSW",
  inherit = BO_base,
  public = list(
    method = 7,
    optimize = function(Budget0=20, Budget=500, num_ref_x=NULL,
                        N0=1000, Na=10, maxevals=100,
                        PN0=4000, PNa=40, Pmaxevals=200, ...){
     
        # define a function that suggests the next x to evaluate
        get_next_x = function(){
          Xr = Build_ref_X(self$GP, T, num_ref_x)
          newx = MCMC_PWKG_grad(list(self$GP), Xr=Xr,
                               N0=N0, Na=Na, maxevals=maxevals, 
                               PN0=PN0, PNa=PNa, Pmaxevals=Pmaxevals)
          return(newx)
        }
     
       # call the optimizer with the suggestion fuction and the kernel
       self$base_optimize(Budget0, Budget, get_next_x, 2)
     }
  )
)


RANDOM_SEARCH = R6Class("RandS", inherit = OptimizerBase,
  public = list(
    x = NULL,
    y = NULL,
    BOseed = 1,
    optimize = function(Budget=500, Ns = 10, BOseed=1,...){
      
      set.seed(BOseed)
      dims = ncol(self$ran)
      ddx = self$ran[2,] - self$ran[1,]
      X = sapply(1:dims, function(d)runif(Budget, self$ran[1,d], self$ran[2,d]))
      S = sample(1:Ns, replace = T, size = Budget)
      XS = cbind(X, S)
      Y  = self$TestFun(XS)
      
      RecX_i = sapply(1:Budget, function(n)which.max(Y[1:n]))
      
      uI = unique(RecX_i)
      uY = TestFun(cbind(X[uI,], 0))
      map = sapply(RecX_i, function(ri)which(ri==uI))
      P = uY[map]
      
      self$Cost = data.frame(N=1:Budget, P=P)
      self$RecX = lapply(1:Budget, function(n)X[RecX_i[n],])
      self$x = X
      self$y = Y
    }
  )
)




