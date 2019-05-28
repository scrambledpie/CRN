source('CRN_BO/GP_model.R')

#############################################################################################
#############################################################################################
## UTILITIES


KGCBcpp   = function(a, b){
  
  if(all(abs(b)<0.000001))return(0)
  
  big = abs(b)>0.000001
  if(any(!big)){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
  }
  
  nn=length(a)
  if(nn==1){
    return(0)
  }
  
  O=order(b)
  b=b[O]
  a=a[O]
  
  while (any(b[-1] == b[-nn])){
    Ri = which(b[-1] == b[-nn])[1]
    if(a[Ri] > a[Ri+1]){
      a=a[-(Ri+1)]
      b=b[-(Ri+1)]
    }else if(a[Ri] <= a[Ri+1]) {
      a=a[-Ri]
      b=b[-Ri]
    }
    nn=length(a)
  }
  
  if(nn==1){
    return(0)
  }
  
  O = KGCB_cpp(c(0,a),c(0,b))
  
  Z=c(-Inf,O[-1,2],Inf)
  I=O[,1]
  
  pC=pnorm(Z)
  dC=dnorm(Z)
  
  sum(  a[I]*(pC[-1]-pC[-length(pC)]) - b[I]*(dC[-1]-dC[-length(pC)])  ) - max(a)
}

Build_ref_X = function(GP, rounding=F){
  
  xd   = GP$xd
  XRAN = GP$XRAN
  
  
  
  xd = xd[,1:ncol(XRAN)]
  xd = Check_X(xd, ncol(XRAN), F, "Xr builder1")
  
  Xr = rbind(LHCran(nrow(xd), XRAN),
             xd + rnorm(length(xd)))
  
  
  Xr = sapply(1:GP$dims, function(d){Xr[,d] = pmax(Xr[,d], XRAN[1,d]); pmin(Xr[,d], XRAN[2,d])})
  
  if(rounding){
    Xr = round(Xr)
    
    repeats = rep(F, nrow(Xr))
    for(i in 1:(nrow(Xr)-1)) for (j in (i+1):nrow(Xr))if(all(Xr[i,]==Xr[j,]))repeats[i] = T
    
    while(any(repeats)){
      Xr[repeats,] = round(LHCran(sum(repeats), XRAN))
      repeats = rep(F, nrow(Xr))
      for(i in 1:(nrow(Xr)-1)) for (j in (i+1):nrow(Xr))if(all(Xr[i,]==Xr[j,]))repeats[i] = T
    }
  }
  
  Xr = cbind(Xr, 0)
  
  Check_ref_X(Xr, ncol(XRAN), T, "Xr builder2")
}

Check_ref_X = function(Xr, dims, withseeds, fname){
  
  Xr = Check_X(Xr, dims, withseeds, fname)
  
  NN = nrow(Xr)
  
  if(nrow(Xr)>1){
    repeats = sapply(1:(NN-1), function(i)any(sapply((i+1):NN, function(j)all(Xr[i,]==Xr[j,]))))
    
    repeats = c(repeats, F)
    
    Xr = Xr[!repeats,,drop=F]
    
  }
  
  Xr = Check_X(Xr, dims, withseeds, paste(fname, "reps removed"))
  
  Xr
  
}

KGCBfilter = function(a,b){
  
  if(all(abs(b)<0.000001))return(0)
  
  big = abs(b)>0.000001
  if(any(!big)){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
  }
  
  
  nn=length(a)
  if(nn==1){
    return(0)
  }
  
  O=order(b)
  b=b[O]
  a=a[O]
  #sort (a,b) by ascending b
  
  
  while (any(b[-1] == b[-nn])){
    Ri = which(b[-1] == b[-nn])[1]
    if(a[Ri] > a[Ri+1]){
      a=a[-(Ri+1)]
      b=b[-(Ri+1)]
    }else if(a[Ri] <= a[Ri+1]) {
      a=a[-Ri]
      b=b[-Ri]
    }
    nn=length(a)
  }
  
  if(nn==1){
    return(0)
  }
  
  A     = 1
  Alast = 1
  C=-Inf
  # print(c(N0(),"alast",Alast,"nn", nn))
  tt =1
  while(Alast<nn){
    s1     = Alast 
    sI     = (s1+1):nn
    CI     = -(a[s1]-a[sI])/(b[s1]-b[sI])
    
    bestsI = which.min(CI)
    Alast  = sI[bestsI]
    A      = c(A,Alast)
    C      = c(C,CI[bestsI])
    # print(c(tt, N0(),"Alast", Alast, "nn", nn, "bestS", bestsI, "as1", a[s1], "asI", a[sI]))
    tt = tt + 1
  }
  C=c(C,Inf)
  
  pC=pnorm(C)
  dC=dnorm(C)
  
  pC_diff = (pC[-1]-pC[-length(pC)])
  
  dC_diff = (dC[-1]-dC[-length(dC)])
  
  grad_a = rep(0, length(a))
  grad_b = rep(0, length(a))
  
  grad_a[A] = pC_diff
  
  sum(  a[A]*pC_diff - b[A]*dC_diff ) - max(a)
  
}

KGCBfilter_grad = function(a,b){
  a_old = a
  
  grad_a = rep(0, length(a))
  grad_a[which.max(a)] = -1
  
  grad_b = rep(0, length(a))
  
  if(all(abs(b)<0.000001))return(0)
  
  big = abs(b)>0.000001
  if(any(!big)>0){
    a   = c( a[big], max(a[!big]) )
    b   = c( b[big], 0)
    new_order = c(which(big), which(a_old==a[length(a)]))
  }else{
    new_order = 1:length(a)
  }
  
  nn=length(a)
  
  if(nn==1){
    return(0)
  }
  
  
  O=order(b)
  b=b[O]
  a=a[O]
  new_order = new_order[O]
  #sort (a,b) by ascending b
  
  while (any(b[-1] == b[-nn])){
    Ri = which(b[-1] == b[-nn])[1]
    if(a[Ri] > a[Ri+1]){
      a=a[-(Ri+1)]
      b=b[-(Ri+1)]
      new_order = new_order[-(Ri+1)]
    }else if(a[Ri] <= a[Ri+1]) {
      a=a[-Ri]
      b=b[-Ri]
      new_order=new_order[-Ri]
    }
    nn=length(a)
  }
  
  
  if(nn==1){
    return(0)
  }
  A     = 1
  Alast = 1
  C=-Inf
  # print(c(N0(),"alast",Alast,"nn", nn))
  tt =1
  while(Alast<nn){
    s1     = Alast 
    sI     = (s1+1):nn
    CI     = -(a[s1]-a[sI])/(b[s1]-b[sI])
    
    bestsI = which.min(CI)
    Alast  = sI[bestsI]
    A      = c(A, Alast)
    C      = c(C, CI[bestsI])
    # print(c(tt, N0(),"Alast", Alast, "nn", nn, "bestS", bestsI, "as1", a[s1], "asI", a[sI]))
    tt = tt + 1
  }
  C=c(C,Inf)
  
  pC=pnorm(C)
  dC=dnorm(C)
  
  a_diff = c(0, a[A]) - c(a[A], 0)
  b_diff = c(0, b[A]) - c(b[A], 0)
  
  # sum(  a[A]*pC_diff - b[A]*dC_diff ) - max(a)
  KG_result = sum(  a_diff*pC - b_diff*dC ) - max(a)
  
  # g_1 = pC[2] #+ 2*C[1]*dC[1
  # print(A)
  # g_2 = pC[3] - pC[2] -1 # + 2* (C[2]*dC[2] - C[1]*dC[1])
  
  grad_a[new_order[A]] = grad_a[new_order[A]] + pC[-1] - pC[-length(pC)]
  
  dC_diff = dC[-1] - dC[-length(dC)]
  C[1] = 0
  C[length(C)] = 0
  
  grad_b[new_order[A]] = -dC[-1] + dC[-length(dC)]
  # browser()
  
  return(list(KG=KG_result, dmu=grad_a, dsig=grad_b))
  
}

UniformDesign_X = function(N0, ran, TestFun, double=0.5, rounding=F, Ns=NULL){
  
  # Sample 'double' of the x points twice, so 1 / (1+dbl) of budget is unique x
  Nx = ceiling(N0/(1+double))
  
  XX = LHCran(Nx, ran)
  
  dims = ncol(ran)
  
  XX = XX[c(sample(1:Nx), sample(1:Nx,N0-Nx)), ,drop=F]
  
  if(rounding)XX = round(XX)
  
  # browser()
  # OO = order(XX[,1])
  # XX = XX[OO,,drop=F]
  
  if(is.null(Ns)) SS = 1:N0 else SS = rep(1:Ns, len=N0)
  
  XX = cbind(XX, SS)
  
  for(i in 1:(N0-1)){
    repped = sapply((i+1):N0, function(j)all(XX[i,]==XX[j,]))
    while(any(repped)){
      XX[i,dims+1] = sample(which(1:Ns!=XX[i,dims+1]), size=1)
      repped = sapply((i+1):N0, function(j)all(XX[i,]==XX[j,]))
    }
  }
  
  
  XX
}

UniformDesign = function(N0, ran, TestFun, double=0.5, rounding=F, Ns=NULL){
  
  # Sample 'double' of the x points twice, so 1 / (1+dbl) of budget is unique x
  Nx = floor(N0/(1+double))
  
  XX = LHCran(Nx, ran)
  
  dims = ncol(ran)
  
  XX = XX[c(1:Nx, sample(1:Nx,N0-Nx)), ,drop=F]
  
  # browser()
  OO = order(XX[,1])
  XX = XX[OO,,drop=F]
  
  if(is.null(Ns)) SS = 1:N0 else SS = rep(1:Ns, len=N0)
  
  XX = cbind(XX, SS)
  
  if(rounding)XX = round(XX)
  
  YY  = TestFun(XX)
  
  if(is.null(Ns)){
    list(x=XX[,1:dims, drop=F], y=YY)
  }else{
    list(x=XX, y=YY)
  }
  
}

CRNUniformDesign = function(Nx, Ns, XRAN, TestFun, rounding=F){
  
  XX = LHCran(Nx, XRAN)
  XX   = XX[rep(1:Nx, Ns),,drop=F]
  
  XXs  = cbind(XX, rep(1:Ns, each=Nx))
  if(rounding)XXs = round(XXs)
  
  YY   = TestFun(XXs)
  list(x=XXs, y=YY)
}

CRNUniformDesign_ATO = function(N0, Ns, XRAN, TestFun, rounding=F){
  
  XX = LHCran(N0, XRAN)
  
  XXs  = cbind(XX, rep(1:Ns, len=N0))
  
  if(rounding)XXs = round(XXs)
  YY   = TestFun(XXs)
  list(x=XXs, y=YY)
}

New_CRNUniformDesign = function(N0, XRAN, TestFun, split=0.5, rounding=F){
  
  N_s = 3
  
  N_x_d = round(N0*split/N_s)
  N_x = N0 - N_x_d*N_s
  
  # Get a set of X points
  XX = LHCran(N_x, XRAN)
  
  # cherry pick some for all seed sampling
  d_i   = sample(N_x, N_x_d)
  XX_d  = XX[d_i, ,drop=F]
  XX_d  = XX_d[rep(1:N_x_d, N_s), ,drop=F]
  XX_d  = cbind(XX_d, rep(1:N_s, each=N_x_d))
  
  # get the remaining X for individual seeds
  s_i   = (1:N_x)[! d_i %in% (1:N_x)]
  XX_s  = XX[s_i, ,drop=F]
  XX_s  = cbind(XX_s, rep(1:N_s, len=nrow(XX_s)))
  
  
  # append x values and 
  XX   = rbind(XX_d, XX_s)
  
  if(rounding) XX = round(XX)
  
  YY   = TestFun(XX)
  list(x=XX, y=YY)
}



#############################################################################################
#############################################################################################
## ACQUISITION FUNCTIONS


Make_CRNKG_grad = function(CRNGP, Xr){
  
  # params:
  #  CRNGP: GP object from GP_model.R
  #  Xr: a set of candidate reference points
  #
  # returns:
  #  list(KG, dKG): the KG acq fun and its derivitive
  #
  
  
  # Check reference points and remove repeats
  Xr = Check_ref_X(Xr, CRNGP$dims, T, "CRNKG Xr")
  ref_seed = Xr[1,CRNGP$dims+1]
  
  # browser()
  iKr = CRNGP$kernel(Xr, CRNGP$xd)%*%CRNGP$iK
  Mr  = CRNGP$MU(Xr)
  
  D = 1:CRNGP$dims
  tXr = t(Xr[D,])
  
  
  
  CRNKG = function(xs){
    
    # repeats = apply(xs[D]==tXr, 2, all)
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats] }
    
    xs = matrix(xs,1)
    # xs = Check_X(xs, CRNGP$dims, T, "CRNKG xs")
    
    # if x is in the set Xr
    # repeats = apply(Xr, 1, function(xi)all(xi[D]==xs[D]))
    
    
    SDx = sqrt(abs(CRNGP$COV(xs, xs)))[1]
    
    # browser()
    
    P0 = CRNGP$kernel(Xr, xs)
    # Nr = P0 > 0.001
    # 
    SIGT = rep(0, length(P0))
    # 
    # SIGT[Nr] = ( P0[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, xs) ) / SDx
    SIGT = ( P0 - iKr%*%CRNGP$kernel(CRNGP$xd, xs) ) * (1/SDx)
    
    xs0  = matrix(c(xs[D], ref_seed), 1)
    
    SIGT = c(SIGT, CRNGP$COV(xs0, xs)/SDx)
    
    MM   = c(Mr, CRNGP$MU(xs0))
    
    O = KGCBfilter(MM, SIGT)
    
    O
  }
  
  dCRNKG = function(xs){
    # xs = Check_X(xs, CRNGP$dims, T, "CRNKG xs")
    xs = matrix(xs, 1)
    # if x is in the set Xr
    # repeats = apply(Xr, 1, function(xi)all(xi[D]==xs[D]))
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats] }
    
    SDx = sqrt(abs(CRNGP$COV(xs, xs)))[1]
    
    # browser()
    
    P0 = CRNGP$kernel(Xr, xs)
    Nr = P0 > -Inf
    # 
    SIGT = rep(0, length(P0))
    # 
    SIGT[Nr] = ( P0[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, xs) ) / SDx
    # SIGT = ( P0 - iKr%*%CRNGP$kernel(CRNGP$xd, xs) ) * (1/SDx)
    
    xs0  = matrix(c(xs[D], ref_seed), 1)
    
    SIGT = c(SIGT, CRNGP$COV(xs0, xs)/SDx)
    
    MM   = c(Mr, CRNGP$MU(xs0))
    
    O = KGCBfilter_grad(MM, SIGT)
    
    Nr = abs(O$dsig)>0
    Nrr = which(Nr[-length(Nr)])
    
    Nr = which(Nr)
    
    # browser()
    
    dKG_s = (O$dsig[Nr]%*%CRNGP$dSIGT.dx_xr(xs, Xr[Nrr,,drop=F], t(iKr[Nrr,,drop=F]))  )[1:CRNGP$dims]
    
    dKG_m = unlist( CRNGP$DTheta(xs[1:CRNGP$dims])) * O$dmu[length(O$dmu)]
    # dKG_s = (O$dsig%*%CRNGP$dSIGT.dx_xr(xs, Xr, t(iKr)))[1:CRNGP$dims]
    
    list(dKG=dKG_m + dKG_s, KG=O$KG)
  }
  
  return(list(KG=CRNKG, dKG=dCRNKG))
}

MCMC_CRNKG_grad = function(CRNGPs, Xr, check_Seeds=NULL,
                           N0=1000, Na=1, maxevals=50){
  
  # params:
  #  CRNGPs: a list of GP objects from GP_model.R
  #  Xr: a set of candidate reference points
  #  N0: number of random search points for AF optimizer
  #  Na: number of the random starts to continue for grad ascent
  #  maxevals: number of grad ascent steps to take
  #
  # returns:
  #  list(KG, dKG): the KG acq fun and its derivitive
  #
  
  
  dims = CRNGPs[[1]]$dims
  
  if(is.null(check_Seeds)) check_Seeds = max(CRNGPs[[1]]$xd[,dims+1])
  
  cat("optim CRNKG_grad, check seeds:", check_Seeds, "...")
  
  # Check reference points and remove repeats
  Xr = Check_ref_X(Xr, CRNGPs[[1]]$dims, T, "CRNKG Xr")
  
  topKG = -Inf
  topxs = 0
  
  if(length(CRNGPs)==1){
    KG_FUNS = Make_CRNKG_grad(CRNGPs[[1]], Xr)
    KG_i = KG_FUNS$KG
    dKG_i = KG_FUNS$dKG
    KG = function(xs){
      KG_v = KG_i(xs)
      if(KG_v > topKG){
        
        topKG <<- KG_v
        topxs <<- xs
        
        # cat("best x ", signif(topxs,3), "  value ", signif(topKG),"\n")
      }
      KG_v
    }
    
    dKG = function(xs){
      O = dKG_i(xs)
      if(O$KG > topKG){
        topKG <<- O$KG
        topxs <<- xs
      }
      return(O$dKG)
    }
    
  }else{
    KG_FUNS = lapply(CRNGPs, function(CRNGP)Make_CRNKG_grad(CRNGP, Xr) )
    
    KG = function(xs){
      KG_v = mean(sapply(KG_FUNS, function(KG_F)KG_F$KG(xs)))
      if(KG_v > topKG){
        topKG <<- KG_v
        topxs <<- xs
      }
      KG_v
    }
    
    dKG = function(xs){
      KG_v = sapply(KG_FUNS, function(KG_F)KG_F$dKG(xs))
      KG_v = apply(KG_v, 1, mean)
      KG_v
    }
    
  }
  
  
  # KG_i = sapply(check_seeds, function(s)sapply(1:99,function(xi)KG(c(xi,s))))
  # plot(c(0,100), range(KG_i), col="white")
  # BB = sapply(check_seeds, function(s)lines(1:99, KG_i[,s], col=s+1))
  
  for(s in check_Seeds){
    KGs = function(x)KG(c(x,s))
    dKGs = function(x)dKG(c(x,s))
    dead = Optimizer_2(KGs, dKGs, ran=CRNGPs[[1]]$XRAN,
                       N0=N0,
                       Na=Na,
                       maxevals=maxevals)
  }
  
  # cat("best KG ", topxs, "\t", topKG,"\t", max(KG_i), "\n")
  cat("done, ")
  topxs
}




Make_PWKG_grad = function(CRNGP, Xr, ratio=1){
  
  # params:
  #  CRNGP: GP object from GP_model.R
  #  Xr: a set of candidate reference points
  #  ratio: a multiplier for the PWKG to turn it on/off
  #
  # returns:
  #  KG: the KG acq fun
  #  dKG: grad of KG
  #  PWKG: pairwise KG
  #  dPWKG: grad of pairwise-KG
  #  best_x: function that returns best observed x
  #
  
  
  Xr = Check_ref_X(Xr, CRNGP$dims, T, "PWKG Xr")
  new_seed = max(CRNGP$xd[,CRNGP$dims+1]) + 1
  
  iKr = CRNGP$kernel(Xr, CRNGP$xd)%*%CRNGP$iK
  Mr  = CRNGP$MU(Xr)
  
  D = 1:CRNGP$dims
  
  tXr = t(Xr[D,])
  
  KGFUNS = Make_CRNKG_grad(CRNGP, Xr)
  
  topKG = -Inf
  topX  = 0
  
  KG = function(x){
    xs = matrix(c(x, new_seed), 1)
    # xs = Check_X(xs, CRNGP$dims, T, "PWKG KG xs")
    O = KGFUNS$KG(xs)
    
    if(O>topKG){
      topKG <<- O
      topX <<-  xs
    }
    
    O
  }
  
  dKG = function(x){
    xs = matrix(c(x, new_seed), 1)
    # xs = Check_X(xs, CRNGP$dims, T, "PWKG KG xs")
    
    O = KGFUNS$dKG(xs)
    
    if(O$KG>topKG){
      topKG <<- O$KG
      topX <<-  xs
    }
    
    return(O$dKG)
  }
  
  PWKG = function(xs1xs2){
    
    x1 = matrix( c(xs1xs2[D], new_seed) , 1)
    x2 = matrix( c(xs1xs2[CRNGP$dims + D], new_seed), 1)
    
    # repeats = apply(x1[D]==tXr, 2, all) | apply(x2[D]==tXr, 2, all)
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats] }
    
    # x1 = Check_X(x1, CRNGP$dims, T, "PWKG PW x1")
    # x2 = Check_X(x2, CRNGP$dims, T, "PWKG PW x2")
    
    if (all(x1==x2)) return(0)
    
    
    SDx = sqrt(abs(abs(CRNGP$COV(x1, x1)) + abs(CRNGP$COV(x2, x2)) - 2*CRNGP$COV(x1, x2)))[1]
    
    P1 = CRNGP$kernel(Xr, x1)
    P2 = CRNGP$kernel(Xr, x2)
    Nr = P1 > 0.001 | P2 > 0.001
    
    SIGT = rep(0, length(P1))
    SIGT[Nr] =            P1[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x1) 
    SIGT[Nr] = SIGT[Nr] - P2[Nr] + iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x2)
    
    xs01  = matrix(c(x1[D], 0), 1)
    xs02  = matrix(c(x2[D], 0), 1)
    
    SIGTx1 = CRNGP$COV(xs01, x1)[1] - CRNGP$COV(xs01, x2)[1]
    SIGTx2 = CRNGP$COV(xs02, x1)[1] - CRNGP$COV(xs02, x2)[1]
    
    SIGT = c(SIGT, SIGTx1, SIGTx2) *(1 / SDx)
    
    MM   = c(Mr, CRNGP$MU(xs01), CRNGP$MU(xs02))
    
    # browser()
    
    O = 0.5 * KGCBfilter(MM, SIGT) * ratio
    
    if(O>topKG){
      topKG <<- O
      topX <<-  rbind(x1, x2)
    }
    
    O
    
  }
  
  dPWKG = function(xs1xs2){
    
    x1 = matrix( c(xs1xs2[D], new_seed), 1)
    x2 = matrix( c(xs1xs2[CRNGP$dims + D], new_seed), 1)
    
    # x1 = Check_X(x1, CRNGP$dims, T, "PWKG PW x1")
    # x2 = Check_X(x2, CRNGP$dims, T, "PWKG PW x2")
    
    if (all(x1==x2)) return(rep(0, CRNGP$dims*2))
    
    # if(xs[1]==75)browser()
    # if x is in the set Xr
    # repeats = apply(Xr, 1, function(xi)all(xi[D]==x1[D])) | apply(Xr, 1, function(xi)all(xi[D]==x2[D]))
    # if (any(repeats)){Xr = Xr[!repeats,]; iKr = iKr[!repeats,]; Mr = Mr[!repeats]}
    # browser()
    VARx = abs(abs(CRNGP$COV(x1, x1)) + abs(CRNGP$COV(x2, x2)) - 2*CRNGP$COV(x1, x2) )[1]
    SDx  = sqrt(VARx)
    
    P1 = CRNGP$kernel(Xr, x1)
    P2 = CRNGP$kernel(Xr, x2)
    Nr = P1 > 0.001 | P2 > 0.001
    
    COV1     = rep(0, length(P1))
    COV1[Nr] = P1[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x1)
    
    COV2     = rep(0, length(P1))
    COV2[Nr] = P2[Nr] - iKr[Nr,]%*%CRNGP$kernel(CRNGP$xd, x2)
    
    xs01  = matrix( c(x1[D], 0), 1)
    xs02  = matrix( c(x2[D], 0), 1)
    
    COV1 = c( COV1, CRNGP$COV(x1, xs01)[1], CRNGP$COV(x1, xs02)[1])
    COV2 = c( COV2, CRNGP$COV(x2, xs01)[1], CRNGP$COV(x2, xs02)[1])
    
    COV_dif = COV1 - COV2
    
    SIGT = COV_dif * (1 / SDx)
    
    MM   = c(Mr, CRNGP$MU(xs01), CRNGP$MU(xs02))
    
    O = lapply(KGCBfilter_grad(MM, SIGT), function(f)ratio*0.5*f)
    
    if(O$KG>topKG){
      topKG <<- O$KG
      topX <<-  rbind(x1,x2)
    }
    
    # browser()
    Nr = Nr & abs(O$dsig[1:(length(O$dsig)-2)])>0
    
    # browser()
    
    ul = function(L)sapply(L,I)
    
    # First dKG/dmu * dmu/dx12
    dkg.dmu = O$dmu[length(O$dmu) + rep(c(-1, 0),each=CRNGP$dims)]
    dMM = ul( c(CRNGP$DTheta(x1[D]), CRNGP$DTheta(x2[D])) )
    
    # Second dKG/dsig * dsig/dx12
    # sigt(x1, x2; Xr,x1,x2) = COV(x1, (Xr,x10,x20)) - COV(x2, (Xr,x10,x20)) / sqrt(SDx)
    # A=numerator dA/dx1 then d/dx2
    XrN = Xr[Nr,,drop=F]
    iKN = t(iKr[Nr,,drop=F])
    
    dCOV1 = matrix(0, nrow(Xr), CRNGP$dims)
    if(sum(Nr)>0) dCOV1[Nr,] = ul(CRNGP$dCOV.dx1(x1, XrN, iKN))
    dCOV1 = rbind(dCOV1, 
                  ul(CRNGP$dCOVxx0.dx(x1)) - ul(CRNGP$dCOV.dx1(xs01, x2)), 
                  ul(CRNGP$dCOV.dx1(x1, xs02)))
    
    
    dCOV2 = matrix(0, nrow(Xr), CRNGP$dims)
    if(sum(Nr)>0) dCOV2[Nr,] = ul(CRNGP$dCOV.dx1(x2, XrN, iKN))
    dCOV2 = rbind(dCOV2, 
                  ul(CRNGP$dCOV.dx1(x2, xs01)), 
                  ul(CRNGP$dCOVxx0.dx(x2)) - ul(CRNGP$dCOV.dx1(xs02, x1)))
    
    dCOV = cbind(dCOV1, -dCOV2)
    
    # B=denominator^2 d/dx COV(x1,x1) + COV(x2,x2) - 2 COV(x1,x2)
    
    dSD1 = ul(CRNGP$dCOVxx.dx(x1)) - 2*ul(CRNGP$dCOV.dx1(x1,x2))
    dSD2 = ul(CRNGP$dCOVxx.dx(x2)) - 2*ul(CRNGP$dCOV.dx1(x2,x1))
    
    dSD = c(dSD1, dSD2)
    
    # dsig/dx12 = (dA - dB A /2B) / sqrt(B)
    
    dSIGT = (dCOV - outer(COV_dif, dSD)*(0.5/VARx)) *(1/SDx)
    
    dPWKG = dkg.dmu*dMM + as.numeric(O$dsig%*%dSIGT)
    
    dPWKG
    
  }
  
  best_x = function(){
    topX
  }
  
  return(list(KG=KG, dKG=dKG, PWKG=PWKG, dPWKG=dPWKG, best_x=best_x))
}

MCMC_PWKG_grad = function(CRNGPs, Xr, ratio=1, 
                          N0=1000, Na=1, maxevals=50, 
                          PN0=4000, PNa=4, Pmaxevals=100){
  
  # params:
  #  CRNGPs: a list of GP objects from GP_model.R
  #  Xr: a set of candidate reference points
  #  ratio: multiplier for PW-KG value
  #  N0: number of random search points for KG optimizer
  #  Na: number of the random starts to continue for grad ascent
  #  maxevals: number of grad ascent steps to take
  #  PN0: number of random search points for PW-KG optimizer
  #  PNa: number of the random starts to continue for grad ascent
  #  Pmaxevals: number of grad ascent steps to take  
  #
  # returns:
  #  list(KG, dKG): the KG acq fun and its derivitive
  #
  
  
  stopifnot(
    typeof(CRNGPs)=="list",
    ncol(Xr)==CRNGPs[[1]]$dims+1
  )
  
  cat("optim PWKG_grad with pairs, ...")
  topKG = -Inf
  topxs = 0
  topPW = FALSE
  
  
  new_seed = max(CRNGPs[[1]]$xd[,CRNGPs[[1]]$dims+1]) + 1
  
  if(length(CRNGPs)==1){
    
    FUNS = Make_PWKG_grad(CRNGPs[[1]], Xr, ratio=ratio)
    
    KG = function(xs){
      O = FUNS$KG(xs)
      if (O>topKG){
        topKG <<- O
        topxs <<- c(xs, new_seed)
        topPW <<- F
        # cat("best single ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dKG = function(xs){
      FUNS$dKG(xs)
    }
    
    
    PWKG = function(xs1xs2){
      O = FUNS$PWKG(xs1xs2)
      if (O>topKG){
        topKG <<- O
        topxs <<- cbind(matrix(xs1xs2, 2, byrow = T), new_seed)
        topPW <<- T
        # cat("best pair ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dPWKG = function(xs1xs2){
      FUNS$dPWKG(xs1xs2)
    }
    
  }else{
    
    FUNS = lapply(CRNGPs, function(CRNGP)Make_PWKG_grad(CRNGP, Xr) )
    
    KG = function(xs){
      O = mean(sapply(FUNS, function(FF)FF$KG(xs)))
      if (O>topKG){
        topKG <<- O
        topxs <<- c(xs, new_seed)
        topPW <<- F
        # cat("best single ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dKG = function(xs){
      dKG_v = sapply(FUNS, function(FF)FF$dKG(xs))
      O = apply(dKG_v, 1, mean)
      O
    }
    
    PWKG = function(xs1xs2){
      O = mean(sapply(FUNS, function(FF)FF$PWKG(xs1xs2)))
      if (O>topKG){
        topKG <<- O
        topxs <<- cbind(matrix(xs1xs2, 2, byrow = T), new_seed)
        topPW <<- T
        # cat("best pair ", topxs, " value ", topKG, "\n")
      }
      O
    }
    dPWKG = function(xs){
      dKG_v = sapply(FUNS, function(FF)FF$dPWKG(xs))
      O = apply(dKG_v, 1, mean)
      O
    }
    
  }
  # par(mfrow=c(1,1))
  # XX = 1:100
  # XXXX = cbind(XX, rep(XX, each=100))
  # KG_eval = sapply(XX, function(xi)KG(c(xi, new_seed)))
  # PWKG_eval = matrix(apply(XXXX, 1, PWKG), 100, 100)
  # contour(XX, XX, PWKG_eval)
  # PWKG_eval = apply(PWKG_eval, 1, max)
  # # plot(range(XX), range(c(KG_eval, PWKG_eval)), col="white")
  # scale = 100 / max(c(KG_eval, PWKG_eval))
  # lines(XX, scale*KG_eval, col="red")
  # lines(XX, scale*PWKG_eval, col="blue", lwd=2)
  
  # browser()
  
  # A = EI_optimizer(KG, CRNGPs[[1]]$XRAN)
  
  # A = EI_optimizer(PWKG, cbind(CRNGPs[[1]]$XRAN, CRNGPs[[1]]$XRAN))
  
  A = Optimizer_2(KG, dKG, CRNGPs[[1]]$XRAN,
                  N0=N0,
                  Na=Na,
                  maxevals=maxevals)
  
  A = Optimizer_2(PWKG, dPWKG, cbind(CRNGPs[[1]]$XRAN, CRNGPs[[1]]$XRAN),
                  N0 = PN0,
                  Na = PNa,
                  maxevals = Pmaxevals)
  
  
  topxs = FUNS$best_x()
  
  if (length(topxs)>CRNGPs[[1]]$dims+1){
    cat("best is pair, ")
  }else{
    cat("best is singleton, ")
  }
  
  return(topxs)
  
}




cat(" Loaded Acquitision_funs.R\n\n")


