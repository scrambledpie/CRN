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


Kernel_regression = function(X_test, X_train, Y_train, r, lb, ub){
    # Makes predictions of output at X_test using top-hat kernel
    #
    # ARGS:
    #  X_test: (n_test, dims) matrix
    #  X_train: (n_train, dims) matrix
    #  Y_train: (n_train) vector
    #  r: radius of kernel in units of (ub - lb)
    #  lb: (dims) lower bounds on X
    #  ub: (dims) upper bounds on X
    #
    # RETURNS:
    #  pred: (n_test) vector

    dims = ncol(X_train)
    n_test = nrow(X_test)
    n_train = nrow(X_train)

    # (dims) float
    invlx_sq = 1/(ub - lb)^2

    # (dims, n_test, n_train)
    diffs_sq = lapply(1:dims, function(i)outer(X_test[ ,i], X_train[ ,i], "-")^2 * invlx_sq[i])

    # (n_test, n_train) float
    x_sq = Reduce("+", diffs_sq)

    # scalar
    r_sq = r * r

    # (n_test, n_train) binary -> float
    mask_mat = 1.0 * (x_sq < r_sq)

    # (n_test, n_train) float
    Y_train_mat = matrix(Y_train, n_test, n_train, byrow=T)

    # (n_test)
    pred = apply( mask_mat * Y_train_mat, 1, sum) / apply( mask_mat, 1, sum)

    return(pred)
}

Root_Finding = function(fun, lo, hi, iters=10){
  # Applies 10 iterations of bisection
  # ARGS:
  #   fun: callable scalar input function
  #   lo: starting lower bound
  #   hi: starting upper bound
  #   iters: number of bisections to perform

  # RETURNS:
  #   lo: new lo bound
  #   hi: new hi bound
  #   f_lo: function at lower bound
  #   f_hi: function at upper bound

  f_lo = fun(lo)
  f_hi = fun(hi)

  stopifnot(
    lo < hi,
    sign(f_lo) != sign(f_hi)
  )

  # plot(c(lo, hi), c(f_lo, f_hi))
  # lines(c(lo, hi), c(0, 0))

  t = 0
  while(t < iters){

    mid = (lo + hi) * 0.5
    f_mid = fun(mid)

    # points(mid, f_mid, col="green", pch=19)

    
    if(f_mid==0){ 
      # found the root!
      t = iters + 1
      lo = mid
      hi = mid
      f_lo = 0
      f_hi = 0

    }else if(sign(f_lo) == sign(f_mid)){
      # if the middle is the same as lo, update lo
      lo = mid
      f_lo = f_mid

    }else{
      # if the middle is the same as hi, update hi
      hi = mid
      f_hi = f_mid

    }

    t = t+1
  }

  return(list(lo=lo, hi=hi, f_lo=f_lo, f_hi = f_hi))
}

if(1==11){

  fun = function(x) -x ^3 - 2

  lo = -2
  hi = 2

  png("mygraphic.png")

  output = Root_Finding(fun, lo, hi)

  X = seq(lo, hi, len=100)
  F = fun(X)

  lines(X, F)

  dev.off()
  browseURL("mygraphic.png")

}

pointwise_ball_regression = function(X_test, X_train, Y_train, r_series, lb, ub){
    # Makes predictions of output at X_test using top-hat kernel around each point.
    # Each point has a unique radius/kernel width.
    #
    # ARGS:
    #  X_test: (n_test, dims) matrix
    #  X_train: (n_train, dims) matrix
    #  Y_train: (n_train) vector
    #  r_series: (>nx) vector of radii of kernel in units of (ub - lb)
    #  lb: (dims) lower bounds on X
    #  ub: (dims) upper bounds on X
    #
    # RETURNS:
    #  pred: (n_test) vector

    dims = ncol(X_train)
    n_test = nrow(X_test)
    n_train = nrow(X_train)

    stopifnot(
      n_train==length(Y_train),
      n_train<=length(r_series),
      dims==length(lb),
      dims==length(ub)
    )


    # (dims) float
    invlx_sq = 1/(ub - lb)^2

    # (dims, n_test, n_train)
    diffs_sq = lapply(1:dims, function(i)outer(X_test[ ,i], X_train[ ,i], "-")^2 * invlx_sq[i])

    # (n_test, n_train) float
    x_sq_mat = Reduce("+", diffs_sq)

    # broadcast upto (n_test, n_train)
    r_series = r_series[1:n_train]
    r_sq = r_series * r_series
    r_sq_mat = matrix(r_sq, ncol=n_train, nrow=n_test, byrow=TRUE)

    # (n_test, n_train) binary -> float
    mask_mat = 1.0 * (x_sq_mat <= r_sq_mat)

    # (n_test, n_train) float
    Y_train_mat = matrix(Y_train, n_test, n_train, byrow=TRUE)

    ball_counts = apply(mask_mat, 1, sum)

    # (n_test)
    pred = apply( mask_mat * Y_train_mat, 1, sum) / ball_counts

    pred[ball_counts==0] = 0

    

    return(pred)
}


Make_random_point = function(centre, lb, ub){
  # Picks a random point around using hte improved hit and run
  # method. Pick a random trajectory passing through centre and
  # a uniform point along the trajectory.
  # ARGS:
  #   centre: (dims) vector in the search space
  #   lb: (dims) vector of lower bounds
  #   ub: (dims) vector of upper bounds
  
  # RETURNS:
  #   x: (dims) vector

  dims = length(centre)

  stopifnot(
    length(lb)==dims,
    length(ub)==dims,
    all(lb) <= ub
  )

  # Make a random direction by taking the direction of an isotropic Gaussian
  dirn = rnorm(dims)
  dirn = dirn / sqrt(sum(dirn^2))
  dirn = dirn * (ub - lb)


  plot(c(0, 0, 1, 1, 0), c(0, 1, 1, 0, 0), type='l')

  # make trajectory function
  traj = function(lambda)centre + lambda * dirn

  x_lo = traj(-10)
  x_hi = traj(10)
  X_t = cbind(x_lo, x_hi)

  lines(X_t[1,], X_t[2,], lty=2)

  # find feasible range of lambda
  within_bounds = function(lambda){

    x = traj(lambda)

    outcome = -0.5 + all(x >= lb) & all(x <= ub)

    points(x[1], x[2], col= 0.5 + outcome)


    cat("x ", x, "bound ", outcome, "\n")
    return(outcome)

  } 

  lambda_min = Root_Finding(within_bounds, -1, 0)$lo

  lambda_max = Root_Finding(within_bounds, 0, 1)$hi

  x = runif(1, lambda_min, lambda_max)

  points(x[1], x[2], pch=19, col='red')

  # rejection sample until x is inside the space
  while(within_bounds(x) < 0){
    x = runif(1, lambda_min, lambda_max)
    points(x[1], x[2], pch=19, col='red')
  }

  print("made a point")
  return(x)

}

if (1==1){
  # TESTING IMPROVED HIT AND RUN SAMPLER
  png("mygraphic.png")

  lb = c(0, 0)
  ub = c(1, 1)
  centre = c(0.9, 0.9)

  a = Make_random_point(centre, lb, ub)

  dev.off()
  browseURL("mygraphic.png")
}


if(1==11){
  # TESTING OF KERNEL REGRESSION

  X = matrix(1:10, ncol=1)
  Y = rnorm(10)

  X_test = matrix(seq(0, 10, len=100), ncol=1)
  # Y_pred = Kernel_regression(X_test, X, Y, r=0.2, lb=0, ub = 10)

  r_series = rep(2, 10) * ( 0.3 ** (0:9) )
  Y_pred = pointwise_ball_regression(X_test, X, Y, r_series, lb=0, ub = 10)



  png("mygraphic.png")

  plot(X, Y)
  lines(X_test, Y_pred)

  dev.off()
  browseURL("mygraphic.png")
}


######################################################################
# A generic BO base class to be inherited by all BO methods
SOSA = R6Class("SOSA",inherit = OptimizerBase,
  public=list(
    Timing = list(),
    myID   = NULL,
    BOseed = NULL,
    method = NULL,
    Ns0    = NULL,
    r0     = NULL,
    gamma  = NULL,
    beta   = NULL,
    s      = NULL,
    initialize = function(testfun, 
                          ran, 
                          BOseed=NULL, 
                          myID=NULL, 
                          Ns0=5,
                          rounding=F,
                          r0=100,
                          gamma=0.9,
                          beta=0.9,
                          s=0.9
                          ){
      # ARGS:
      #   testfun: callable, takes x and seed, returns scalar
      #   ran: (2, dims) matrix, lower and upper bounds
      #   BOseed: fixed RNG seed
      #   myID: number to use as savefile
      #   Ns0: number of starting points
      #   rounding: bool, stick to integers
      #   r0: initial radius of ball
      #   gamma: order of local sample density
      #   beta: radius decay exponent
      #   s: recomendation time sequence shrinkage

      super$initialize(testfun, ran, rounding)
      self$Ns0     = Ns0
      self$BOseed  = BOseed
      self$myID    = myID
      set.seed(BOseed)

      self$r0 = r0
      self$gamma = gamma
      self$beta = beta
      self$s = s},

    get_RecX = function(){

      if(length(self$predictions)==length(self$Y)){

      }

    },
#
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
#
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