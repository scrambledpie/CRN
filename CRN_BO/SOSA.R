library(R6)


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
  # TESTING ROOT FINDER
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

shrinking_ball_regression = function(X_test, X_train, Y_train, r_series, lb, ub){
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

    stopifnot(
      length(dim(X_train))==2
    )

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

    # model is not defined for x outside of the balls!
    stopifnot(all(ball_counts>0))

    # (n_test)
    pred = apply( mask_mat * Y_train_mat, 1, sum) / ball_counts


    return(pred)
}

if(1==11){
  # TESTING OF KERNEL REGRESSION

  X = matrix(1:10, ncol=1)
  Y = rnorm(10)

  X_test = matrix(seq(0, 10, len=100), ncol=1)
  # Y_pred = Kernel_regression(X_test, X, Y, r=0.2, lb=0, ub = 10)

  r_series = rep(2, 10) * ( 0.3 ** (0:9) )
  Y_pred = shrinking_ball_regression(X_test, X, Y, r_series, lb=0, ub = 10)



  png("mygraphic.png")

  plot(X, Y)
  lines(X_test, Y_pred)

  dev.off()
  browseURL("mygraphic.png")
}

improved_hit_and_run = function(centre, lb, ub){
  # Picks a random point using the improved hit and run
  # method. Pick a straight trajectory passing through centre and
  # with a random direction, pick a uniform point along the trajectory.
  # ARGS:
  #   centre: (dims) point in the search space
  #   lb: (dims) vector of lower bounds
  #   ub: (dims) vector of upper bounds
  #
  # RETURNS:
  #   x: (dims) point in search space

  dims = length(centre)

  stopifnot(
    length(lb)==dims,
    length(ub)==dims,
    all(lb <= ub)
  )

  # Make a random direction from an isotropic Gaussian
  dirn = rnorm(dims)
  dirn = dirn * sqrt(dims) / sqrt(sum(dirn^2))
  dirn = dirn * (ub - lb)

  # make trajectory function
  traj = function(lambda1) (centre + (lambda1 * dirn))

  # indicator function for in/out of search space
  within_bounds = function(lambda2){
    stopifnot(length(lambda2)==1)
    xx = traj(lambda2)
    outcome = -0.5 + 1 * (all(xx >= lb) & all(xx <= ub) )
    return(outcome)
  } 

  # bounds of lambda values, trajectory within space.
  lambda_min = Root_Finding(within_bounds, -1, 0)$lo
  lambda_max = Root_Finding(within_bounds, 0, 1)$hi

  lambda_random = runif(1, lambda_min, lambda_max)

  # rejection sample until lambda_random is inside the space.
  while(within_bounds(lambda_random) < 0){
    lambda_random = runif(1, lambda_min, lambda_max)
  }

  # get the corresponding location
  x = traj(lambda_random)

  return(x)

}

if (1==11){
  #############################################
  # TESTING IMPROVED HIT AND RUN SAMPLER
  pdf("mygraphic.pdf")

  lb = c(0, 0)
  ub = c(1, 1)
  centre = runif(2) #c(0.9, 0.9)

  plot(c(-1.5, 2.5), c(-1.5, 2.5), col="white")
  lines(c(0, 0, 1, 1, 0), c(0, 1, 1, 0, 0))

  a_full = centre
  for( i in 1:1000){
    a = improved_hit_and_run(centre, lb, ub)
    points(a[1], a[2], col='red', pch=4)

    a_full = rbind(a_full, a)
  }

  points(centre[1], centre[2], pch=19, col="brown")
  dev.off()
  browseURL("mygraphic.pdf")
  # system("open -a Preview mygraphic.pdf")

  print(max(a_full))
  print(min(a_full))

}




######################################################################
# A generic BO base class to be inherited by all BO methods
SOSA = R6Class("SOSA",
  public=list(
    TestFun    = NULL,
    ran        = NULL,
    lb         = NULL,
    ub         = NULL,
    RecX       = list(),
    Cost       = data.frame(N=numeric(0), P=numeric(0), slowP=numeric(0)),
    rounding   = FALSE,
    Timing = list(),
    myID   = NULL,
    BOseed = NULL,
    method = NULL,
    r0     = NULL,
    gamma  = NULL,
    beta   = NULL,
    s      = NULL,
    r_series = NULL,
    i_n_series = NULL,
    X      = NULL,
    Y      = NULL,
    Y_pred = NULL,

    initialize = function(testfun, 
                          ran,
                          BOseed=NULL, 
                          myID=NULL, 
                          rounding=F,
                          r0=0.05,
                          gamma=0.9,
                          beta=0.009,
                          s=0.9
                          ){
      # ARGS:
      #   testfun: callable, takes x and seed, returns scalar
      #   ran: (2, dims) matrix, lower and upper bounds
      #   BOseed: fixed RNG seed
      #   myID: number to use as savefile
      #   rounding: bool, stick to integers
      #   r0: initial radius of ball
      #   gamma: order of local sample density
      #   beta: radius decay exponent
      #   s: recomendation time sequence shrinkage

      self$TestFun = TestFun
      self$ran = ran
      self$rounding = rounding
      self$lb = self$ran[1, ]
      self$ub = self$ran[2, ]
      self$BOseed  = BOseed
      self$myID    = myID
      set.seed(BOseed)

      self$r0 = r0
      self$beta = beta
      self$r_series = r0 * ((0:2000) ^ -beta)
      self$gamma = gamma
      self$s = s
      self$i_n_series = floor( 0:2000 ^ s )
    },

    get_RecX = function(){
      # best predicted x in 1,...,n, and best predicted x in 1,..,i_n points

      # if necessary, update the predicted function values
      self$update_pred_Y()
      dims = ncol(self$X) - 1

      # best of all observed X
      rec_X_i = which.max(self$Y_pred)
      best_X = self$X[rec_X_i, 1:dims]

      # best of observed x upto time i_n
      n = length(self$Y)
      i_n = self$i_n_series[n]
      slow_Y_pred = self$Y_pred[1:i_n]
      slow_rec_X_i = which.max(slow_Y_pred)
      slow_best_X = self$X[slow_rec_X_i, 1:dims]

      return(list(best_x=best_X, slow_best_x=slow_best_X))
    },

    get_next_x = function(){
      # Acquisition method: hit and rtun centred 
      # on the predicted peak.
      dims = ncol(self$X) - 1

      # if necessary, update the predicted function values.
      self$update_pred_Y()

      # best of all observed X is the centre of sampler.
      rec_X_i = which.max(self$Y_pred)
      best_X = self$X[rec_X_i, 1:dims]

      new_x = improved_hit_and_run(best_X, self$lb, self$ub)

      # add the seed on
      new_x = c(new_x, length(self$Y))

      return(new_x)
    },

    UpdateLogs = function(KG_time=0, eval_time=0, fit_time=0){

      N = length(self$Y)
      
      self$Timing[[N]] = c(KG_time, eval_time, fit_time)
      self$RecX[[N]]   = self$get_RecX()
      
      Rx_best               = matrix(c(self$RecX[[N]]$best_x, 0), nrow=1)
      Rx_slow               = matrix(c(self$RecX[[N]]$slow_best_x, 0), nrow=1)
      NCi                   = nrow(self$Cost) + 1

      # cat("Rx_best: ", Rx_best, ",  Rx_slow: ", Rx_slow, "\n")
      cat("\nTesting")
      self$Cost[NCi,]  = c(N, self$TestFun(Rx_best), self$TestFun(Rx_slow))
    },

    update_pred_Y = function(){

      stopifnot( length(dim(self$X))==2 )

      if(length(self$Y_pred)!= length(self$Y)){
        dims = ncol(self$X) - 1
        X_train = self$X[ , 1:dims, drop=F]
        self$Y_pred = shrinking_ball_regression(X_train, X_train, self$Y, self$r_series, self$lb, self$ub)
      }
    },

    optimize = function(Budget=500, ...){
      
      
      cat("Starting Optimization")
      t0       = proc.time()[3]
      

      # get first initilialization points
      cat("\nSelecting point 1 \n")
      X_init = ( self$ran[1, ] + self$ran[2, ] ) * 0.5
      if(self$rounding) X_init = round(X_init)
      X_init = matrix(c(X_init, 1), nrow=1) # append the seed
      KG_time = proc.time()[3] - t0
      
      # get objective function values
      Y_init    = self$TestFun(X_init)
      eval_time = proc.time()[3] - t0 - KG_time
      
      # fit the model
      self$X = X_init
      self$Y = Y_init
      self$update_pred_Y()
      fit_time = proc.time()[3] - t0 - KG_time - eval_time
      
      # synchronise status
      self$UpdateLogs(KG_time, eval_time, fit_time)
      
      # Now do the sequential bit!
      while (length(self$Y)<Budget){

        N =  length(self$Y)+1
        
        cat("\n\nSelecting point: ",N, "\n")
        t0 = proc.time()[3]
        
        # Get the new x to evaluate
        new_x = self$get_next_x()
        if(self$rounding) new_x = round(new_x)

        cat("new x is: ", new_x, ", ")
        KG_time = proc.time()[3] - t0
        
        
        # Evaluate the objective function
        new_y = self$TestFun(new_x)
        eval_time = proc.time()[3] - t0 - KG_time
        
        
        # Update the data and either do a full hpar update or just finetune
        self$X      = rbind(self$X, new_x)
        self$Y      = c(self$Y, new_y)
        fit_time    = proc.time()[3] - t0 - eval_time - KG_time
        
        self$UpdateLogs(KG_time, eval_time, fit_time)
        
      }
      
    }
  )
)




if(1==1){
  # RUNNING ON ATO FOR TESTING
  problem  = 2
  BOseed   = 1
  Budget   = 25
  filename = "debug_SOSA" # where to save results

  # Make the TestFunction
  source('TestFuns/TestFuns.R')
  TESTFUNS = c(Build_Xie_ATO_cpp_Testfun,
               Build_Ambulance_Testfun,
               Build_DailyAmbulance_Testfun)

  TestFun = TESTFUNS[[problem]](BOseed, numtestseeds=2000, runlength=1)[[1]]
  ran = attr(TestFun, 'ran') # bounds of the search space.

  print(ran)

  # pick the optimizer from the list
  Optimizer = SOSA$new(TestFun, ran, BOseed, myID=filename, rounding=T)

  Optimizer$optimize(Budget=Budget)
}