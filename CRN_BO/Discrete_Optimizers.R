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
    initialize = function(TestFun, ran, rounding=T, X_domain=NULL){
      self$TestFun = TestFun
      self$ran = ran
      self$rounding = rounding
      self$X_domain = X_domain
    }
  )
)


######################################################################
# A generic BO base class to be inherited by all BO methods
BO_base_disc = R6Class("BO_base_disc",inherit = OptimizerBase,
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
    X_domain=NULL,
    maxY_f = NULL,
    initialize = function(testfun, ran, X_domian, Y_f=NULL, BOseed=NULL,
                          myID=NULL, Ns0=5){
      super$initialize(testfun, ran, rounding=T)
      self$Ns0      = Ns0
      self$BOseed   = BOseed
      self$myID     = myID
      self$X_domain = X_domain
      self$maxY_f      = max(Y_f)
      set.seed(BOseed)
    },
    UpdateLogs = function(KG_time=0, eval_time=0, fit_time=0){
      N = length(self$GP$yd)

      self$GP$Lhood_Prep()
      self$Lhoods[[N]] = self$GP$Lhood_standard(self$GP$HP)

      self$Hpars[[N]]  = self$GP$HP
      self$Ymean[[N]]  = self$GP$ymean
      self$Timing[[N]] = c(KG_time, eval_time, fit_time)
      self$RecX[[N]]   = self$X_domain[which.max(self$GP$MU(cbind(X_domain, 0)))]

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
          self$maxY_f  <<- Input$maxY_f
          
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
                        maxY_f    = self$maxY_f,
                        .Random.seed = .Random.seed
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
    base_optimize = function(Budget0=5, Budget=50, get_next_x=NULL, Ns0=3, hypers=NULL){


      if(is.null(get_next_x))stop("provide a get_next_x() function!")

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
        self$GP$Refresh(learnHpars = 7, Hpars=hypers, ymean=0)
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
        self$GP$Refresh(ymean=0, learnHpars=0)
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

BO_KG_DISC = R6Class("BO_KG_DISC",
  inherit = BO_base_disc,
  public = list(
    method = 2,
    optimize = function(Budget0=5, Budget=50, hypers, ...){

      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        cat(" evaluating KG_DISC, ")
        Xr = cbind(self$X_domain, 0)
        KG = Make_PWKG_grad(self$GP, Xr)$KG
        newx = X_domain[ which.max(sapply(X_domain, KG)) ]
        newx = c(newx, max(self$GP$xd[,ncol(self$GP$xd)])+1)
        return(newx)
      }

      hypers[length(hypers)] = sum(hypers[length(hypers)-(0:2)])
      hypers[length(hypers)-1] = 0
      hypers[length(hypers)-2] = 0

      # call the optimizer with the suggestion fuction and the kernel
      self$base_optimize(Budget0, Budget, get_next_x, hypers=hypers, ...)
    }
  )
)

######################################################################

BO_CRNKG_DISC = R6Class("BO_CRNKG_DISC",
  inherit = BO_base_disc,
  public = list(
    method = 4,
    optimize = function(Budget0=5, Budget=50,...){

      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        cat(" evaluating CRNKG_DISC, ")
        X_domain = self$X_domain
        Xr    = cbind(X_domain, 0)
        CRNKG = Make_CRNKG_grad(self$GP, Xr)$KG

        checkseeds = 1:(max(self$GP$xd[,ncol(self$GP$xd)])+1)
        KGvals = sapply(checkseeds, function(s)sapply(X_domain, function(xi)CRNKG(c(xi, s))))

        newx = X_domain[ which.max(apply(KGvals, 1, max)) ]
        news = which.max(apply(KGvals, 2, max))

        newxs = c(newx, news)

        return(newxs)
      }

      # call the optimizer with the suggestion fuction and the kernel
      self$base_optimize(Budget0, Budget, get_next_x, ...)
    }
  )
)

######################################################################

BO_PWKG_DISC = R6Class("BO_PWKG_DISC",
  inherit = BO_base_disc,
  public = list(
    method = 5,
    optimize = function(Budget0=5, Budget=50, ...){

      # define a function that suggests the next x to evaluate
      get_next_x = function(){
        cat(" evaluating PWKG_DISC, ")
        X_domain = self$X_domain
        Xr       = cbind(X_domain, 0)
        PKG      = Make_PWKG_grad(self$GP, Xr)
        KG_vals  = sapply(X_domain, PKG$KG)

        PW_vals  = sapply(X_domain, function(x1)sapply(X_domain, function(x2)if(x2>=x1) 0 else PKG$PWKG(c(x1,x2))))

        if(max(KG_vals) > max(PW_vals)){
          newx = X_domain[ which.max( KG_vals )]
          newx = c(newx, max(self$GP$xd[,2])+1)
          cat("best is single, ")
        }else{
          newx1 = X_domain[ which.max(apply(PW_vals,1,max)) ]
          newx2 = X_domain[ which.max(apply(PW_vals,2,max)) ]
          newx = cbind(c(newx1, newx2), max(self$GP$xd[,2])+1 )
          cat("best is pair, ")
        }
        return(newx)
      }

     # call the optimizer with the suggestion fuction and the kernel
     self$base_optimize(Budget0, Budget, get_next_x, ...)
   }
  )
)
