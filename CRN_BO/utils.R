#############################################################################################
#  UTILITIES
CPU = system("hostname",intern=T)

# if (grepl("Book", CPU)){
if (F){
  Check_X = function(x, dims, withseeds, fname){
    # Args
    #   x: input vec or matrix
    #   dims: number of continuous dims of x
    #   withseeds: should the x include an integer seed aslast element
    #   fname: printed with the error message
    #
    # Returns
    #   x: correctly formatted as a matrix
    
    dimsN = dims + withseeds
    
    if (dims==1&is.null(dim(x))&!withseeds){ cat(fname, "converting x to matrix"); x = matrix(x, ncol=1)}
    
    if (length(x)==dimsN & length(dim(x))==0){ cat(fname, "converting x to matrix"); x = matrix(x, ncol=dimsN, nrow=1)}
    
    if (length(x)%%dimsN!=0) stop(fname, " wrong length input, got ", length(x))
    
    if (length(dim(x))!=2) stop(fname, "wrong input dim, got ", dim(x))
    
    if (ncol(x)!=dimsN) stop(fname, "wrong # input columns, got ", ncol(x))
    
    if (withseeds & !all(x[,dimsN]%%1==0)) stop(fname, "non-integer seeds not allowed")
    
    if (withseeds & !all(x[,dimsN]>=0)) stop(fname, "negative seeds not allowed")
    
    x
  }
}else{
  Check_X = function(x, dims, withseeds, fname) x
}

LHCdim = function(N0, dims){
  # Args
  #   N0: the number of points to generate
  #   dims: the box for scaling+centring samples
  #   
  # Returns
  #   a randomly generated LHC of N0 points in [0,1]^dims
  
  sapply(1:dims, function(d){  (sample(1:N0) - runif(N0))  })*(1/N0)
}

LHCran = function(N0, ran){
  # Args
  #   N0: the number of points to generate
  #   ran: the box for scaling+centring samples
  #   
  # Returns
  #   a randomly generated LHC of N0 points in box ran dimensions
  
  if(length(ran)==2)ran=matrix(ran,ncol=1)
  
  if (any(ran[1,]>=ran[2,])|nrow(ran)!=2) stop("LHC bounds are not valid: ", paste(bottoms, tops, collapse=", "))
  
  dims = ncol(ran)
  LHC1 = sapply(1:dims,function(d){  (sample(1:N0) - runif(N0))  })*(1/N0)
  
  BB = matrix(ran[1,], N0, length(ran[1,]), byrow=T)
  CC = matrix(ran[2,]-ran[1,], N0, length(ran[1,]), byrow=T)
  
  BB + LHC1*CC
}

Optimizer_2  = function(f1, df1, ran, N0=1000, Na=5, x0=NULL, debugging=F, reltol=1e-8, maxevals=Inf){
  
  
  cat("\n Na:", Na, ", N0: ", N0, ", maxevals: ",maxevals, ",  ")
  # browser()
  # Function that maximizes f1 and returns argmax
  # 
  # Args
  #   f1: callable scalar output continuous function
  #   df1: gradient of f1
  #   ran: upper+lower bounds of all arguments
  #   N0: number of random intial function evaluations
  #   Na: nuber of the random initial points to use for grad ascent
  #   x0: an optional starting point for the optimizer
  # 
  # Returns
  #   bestx: argmax of f1
  
  con = list(reltol = reltol, maxit=75)
  
  bs = ran[1,]
  ts = ran[2,]
  
  stopifnot(!any(bs>=ts),
            !length(bs)!=length(ts),
            !any(is.infinite(bs)),
            !any(is.infinite(ts)),
            !any(is.na(bs)),
            !any(is.na(ts)),
            !(is.null(x0)&N0==0),
            !is.null(bs),
            !is.null(ts)
  )
  
  
  #################################################################################
  # Make random initializer and History
  init = function()runif(length(bs), bs, ts)
  x1   = init()
  
  History = matrix(c(x1, f1(x1)) , 1)
  evals   = 0
  
  #################################################################################
  # objective function wrappers
  ff1 = function(x){
    # cat(evals," ")
    if(evals>=maxevals)stop("reached max evals")
    if(any(x<ran[1,]))stop("out of range, too low! ", paste(x, collapse = ", "))
    if(any(x>ran[2,]))stop("out of range, too high! ", paste(x, collapse = ", "))
    
    # evaluate objective and sanity check
    OO = f1(x)
    
    stopifnot(
      !is.na(OO),
      !is.infinite(OO),
      !is.null(OO)
    )
    
    if(all(x>bs)&all(x<ts)) History <<- rbind(History,c(x,OO))
    
    evals <<- evals+1
    return(-OO)
  }
  
  dff1 = function(x){
    -df1(x)
  }
  
  safe_optim = function(x_start){
    
    evals = 0
    dead = tryCatch({ 
      optim(par = x_start, fn = ff1, dff1, method = "L-BFGS-B", control = con) 
    }, 
    error = function(e)-1e9,
    warning = function(w)-1e9
    )
    
    evals = 0
    dead = tryCatch({ 
      optim(par = x_start, fn = ff1, dff1, method = "CG", control = con) 
    }, 
    error = function(e)-1e9,
    warning = function(w)-1e9
    )
  }
  
  safe_ff1 = function(x){
    tryCatch({ 
      ff1(x)
    }, 
    error = function(e)NaN,
    warning = function(w)NaN
    )
  }
  
  #################################################################################
  # If we have a given starting point, then optimise from that start!
  if (!is.null(x0)){
    if(is.null(dim(x0))) x0 = matrix(x0, 1)
    dead = apply(x0, 1, safe_ff1)
    dead = apply(x0, 1, safe_optim)
  }
  
  #################################################################################
  # If we want to do some random settings
  if (N0>0){
    
    X0      = LHCran(N0, ran)
    evals   = -N0
    Output  = apply(X0, 1, safe_ff1)
    
    Keepers = !is.nan(Output)
    Output  = Output[Keepers]
    X0      = X0[Keepers,,drop=F]
    
    rank    = order(Output, decreasing=T)
    Na      = min(Na, length(Output))
    
    X0      = History[rank[1:Na], , drop=F]
    
    dead    = apply(X0, 1, safe_optim)
    
  }
  
  #################################################################################
  # Get the best and run optimizers one more time
  lHist = History[, length(ts)+1]
  bestx = History[which.max(lHist), 1:length(ts)]
  dead  = safe_optim(bestx)
  
  
  #################################################################################
  # Get the final best and return it
  lHist = History[, length(ts)+1]
  bestx = History[which.max(lHist), 1:length(ts)]
  return( list(xmax = bestx, fmax = max(lHist))  )
  
}

Freezer    = function(f1, df1, ran, fdims, fvals, x0=NULL,...){
  # Function that optimizes f1 with some fixed inputs
  # 
  # Args
  #   f1: callable scalar output continuous function
  #   df1: gradient of f1
  #   ran: lower+upper bounds of all arguments
  #   fdims: list of indices of frozen dims
  #   fvals: values of frozen dimensions
  #   ...: arguments passed to Optimizer()
  # 
  # Returns
  #   bestx: only unfrozen dims argmax of f1
  
  bs = ran[1,]
  ts = ran[2,]
  stopifnot(
    all(bs<ts),
    length(bs)==length(ts),
    all(is.finite(bs)),
    all(is.finite(ts)),
    # !is.na(f1(bs)),
    # !is.na(f1(ts)),
    # !any(is.na(df1(bs))),
    # !any(is.na(df1(ts))),
    !is.null(bs),
    !is.null(ts),
    # !is.null(f1(bs)),
    # !is.null(f1(ts)),
    # length(df1(bs))==length(bs),
    length(fdims)<length(bs),
    length(fvals)==length(fdims),
    # all(bs[fdims]<fvals),
    all(fvals<ts[fdims])
  )
  
  nn      = ncol(ran)
  
  odims = (1:nn)[!1:nn%in%fdims]
  
  ran_o = ran[,odims]
  
  calls = 0
  x1 = 1:nn
  f2  = function(x){
    x1[fdims] = fvals
    x1[odims] = x  
    calls <<- calls+1 
    f1(x1)
  }
  df2 = function(x){
    x1[fdims]=fvals
    x1[odims]=x
    df1(x1)[odims]
  }
  
  if(!is.null(x0)){
    if(length(x0)==ncol(ran)){
      x0 = x0[odims]
    }else if(length(x0)!=length(odims)){
      x0=NULL
    }
  }
  # browser()
  O = Optimizer_2(f2, df2, ran_o, x0=x0, ...)
  
  O
}

EI_optimizer = function(EIfun, ran, Ns=50){
  
  # Function that optimizes EIfun and returns argmax
  # 
  # Args
  #   EIfun: callable scalar output continuous function
  #   ran: lower+upper bounds of all arguments
  #   Ns: number of random optimizer starts 
  # 
  # Returns
  #   bestx: argmax of EIfun
  
  bottoms = ran[1,]
  tops = ran[2,]
  
  if (any(bottoms>=tops)|nrow(ran)!=2) stop("EI optim bounds are not valid: ", paste(bottoms, tops, collapse=", "))
  
  sXA = LHCran(Ns, ran)
  
  topO = -Inf
  topxa = 0
  t=0
  Obj = function(xa){
    if(any(xa<bottoms)|any(xa>tops)) return(0)
    
    out = EIfun(xa)
    if (out > topO){
      topO  <<- out
      topxa <<- xa
    }
    out
  }
  
  if(length(bottoms)>1){
    OO=sapply(1:Ns,function(i)optim(sXA[i,],Obj,method = "Nelder-Mead",control = list(maxit = 75, fnscale = -1)))
  }else{
    
    Xvals = seq(bottoms, tops, len=100)
    best  = which.max(sapply(Xvals[-c(1,100)], Obj)) +1
    IT    = c(Xvals[best-1], Xvals[best+1])
    OO    = optimise(f = Obj,interval = IT,maximum = T)
  }
  
  return( list(x=topxa, fn=topO) )
}

cat(" Loaded utils.R, ")