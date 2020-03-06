library(R6)

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


Kernel_regression = function(X_test, X_train, Y_train, r){
    # Makes predictions of output at X_train

    dims = ncol(X_train)
    n_test = nrow(X_test)

    diffs_sq = lapply(1:dims, function(i)outer(X_test[ ,i], X_train[ ,i], "-")^2)

    x_sq = Reduce("+", diffs_sq)

    r_sq = r * r

    pred = rep(0, n_test)

    for(i in 1:n_test){

        mask_i = x_sq[i, ] < r_sq
        pred[i] = mean( Y_train[mask_i] )

    }

    return(pred)

}


Build_pred_model = function(X_train, Y_train, r_series)



X = matrix(1:10, ncol=1)
Y = rnorm(10)

X_test = matrix(seq(0, 10, len=100), ncol=1)
Y_pred = Kernel_regression(X_test, X, Y, 0.2)



png("mygraphic.png")

plot(X, Y)
lines(X_test, Y_pred)

dev.off()
browseURL("mygraphic.png") 



######################################################################
# A generic BO base class to be inherited by all BO methods
SOSA = R6Class("SOSA",inherit = OptimizerBase,
  public=list(
    Timing = list(),
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

    

    get_RecX = function(){

        n_max = length(self$Y)

        # the prediction model is sum of KNNs at sampled locations
        if(is.null(self$predictions)){

            self$predictions = lapply(n in 1:n_max, function(n){
                model[[n]] = Kernel_regression(self$X[1:n, ], self$X[1:n, ], self$Y[1:n, ], self$r[n])
            }
        }

        # compute the KNN for 
    },

    UpdateLogs = function(KG_time=0, eval_time=0, fit_time=0){
      N = length(self$GP$yd)
      
      self$Timing[[N]] = c(KG_time, eval_time, fit_time)
      
      if(length(self$RecX)>0){
        self$RecX[[N]]   = self$get_RecX()
      }else{
        self$RecX[[N]]   = self$get_RecX()
      }
      
      Rx               = matrix(c(self$RecX[[N]], 0), 1)
      NCi              = nrow(self$Cost)+1
      self$Cost[NCi,]  = c(N, self$TestFun(Rx))
    },

    base_optimize = function(Budget0=20, Budget=500, get_next_x=NULL, learn_kernel=NULL, Ns0=5){
      
      
      if(is.null(get_next_x))stop("provide a get_next_x() function!")
      

        cat("Starting Optimization")
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