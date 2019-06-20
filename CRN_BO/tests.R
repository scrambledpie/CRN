rm(list = ls())

setwd("/Users/pearce/CRN/")

source("CRN_BO/GP_model.R")
source("CRN_BO/Acquisition_funs.R")
library(testit)


tester_gradients = function(f, df, x0, h=NULL, tol=0.01, msg){
  
  # INPUTS
  #  f: callable function, scalar
  #  df: callable function, vector
  #  x0: the point at which to test gradients
  #  h: size of finite difference
  #  tol: threshold for %delta betweren empirical and alanyrtical grads
  #  msg: name of function being tested
  
  if(is.null(h)) h = min( 0.000001, min(abs(x0))/100)
  cat("\nTesting gradients ", msg, "\n")
  
  f0 = f(x0)
  
  ss = 1:length(x0)
  
  emp_df = sapply(ss, function(i)f(x0+h*(i==ss)))
  emp_df = (emp_df-f0)/h
  
  df0 = df(x0)
  
  cat("emp grad:", emp_df, "\n")
  cat("thy grad:", df0, "\n")
  
  too_small = abs(emp_df)<1e-9
  
  emp_df = emp_df[!too_small]
  df0 = df0[!too_small]
  
  assert("gradients dif under tol", max(abs((emp_df-df0)/emp_df))<tol )
  
  cat("gradient ",msg," Test Passed!\n\n")
  
}

################################
# GP model hypers and kernel gradient tests

test_Lhoodgradients = function(){
  
  cat("Testing lhood gradients\n")
  
  # make some toy data
  X_init = matrix(1:10, ncol=1)
  Y_init = sin(X_init)
  ran = matrix(c(0,11),2)
  
  hpars = c(2, 1, 0.01)
  l_hpars = log(hpars)
  h = 0.001
  
  ###########################################################################
  # log search space
  FUNS = LhoodBuilder(X_init, Y_init, F, ran)
  tester_gradients(FUNS$LH, FUNS$DLH, l_hpars, h, msg="log-scale GP lhood")
  
  FUNS = LhoodBuilder(X_init, Y_init, T, ran)
  tester_gradients(FUNS$LH, FUNS$DLH, l_hpars, h, msg="log-scale GP lhood+prior")
  
  lh = FUNS$LH(l_hpars)
  
  ###########################################################################
  # standard search space
  FUNS = LhoodBuilder_standard(X_init, Y_init, F, ran)
  tester_gradients(FUNS$LH, FUNS$DLH, hpars, h, tol=0.1, msg="normal-scale GP lhood")
  
  FUNS = LhoodBuilder_standard(X_init, Y_init, T, ran)
  tester_gradients(FUNS$LH, FUNS$DLH, hpars, h, tol=0.1, msg="normal-scale GP lhood+prior")
  
  lh2 = FUNS$LH(hpars)
  
  ###########################################################################
  
  assert("lhoods differ by more than 1%", abs(lh-lh2)/max(lh,lh2)<0.01 )
  
  
  cat("ALL CRN gradient Tests Passed!\n")
}

test_CRNLhood_grads = function(){
  
  # make some toy data
  set.seed(1)
  seeds = sample(1:10, replace = T, size = 100)
  
  X_init = cbind(rep(1:10, 10), rep(1:10, each=10), seeds)
  Y_init = 20*sin(sqrt(X_init[,1]^2 + X_init[,2]^2))
  ran = matrix(c(0, 11), 2, 2)
  
  hpars = c(1.25, 1.25, 10, 
            1.1, 1.1, 0.1,
            0.1, 0.1)
  
  GP1 = CRNLHood$new(X_init, Y_init, ran, copula = F)
  GP1$Lhood_Prep()
  
  cat("\n")
  
  tester_gradients(f   = GP1$Lhood,
                 df  = GP1$DLhood,
                 x0  = log(hpars),
                 h   = 0.0001,
                 msg = "log scale CRN gradient")
  
  cat("\n")
  tester_gradients(f   = GP1$Lhood_standard,
                 df  = GP1$DLhood_standard,
                 x0  = hpars,
                 h   = 0.001,
                 msg = "CRN gradient")
  
  cat("\n")
  tester_gradients(f   = GP1$LHyperPrior,
                   df  = GP1$DLHyperPrior,
                   x0  = hpars,
                   h   = 0.001,
                   msg = "CRN prior gradient")
  
}

test_CRN_kernel_grad = function(){
  
  # make some toy data
  set.seed(1)
  seeds = sample(1:10, replace = T, size = 100)
  
  X_init = cbind(rep(1:10, 10), rep(1:10, each=10), seeds)
  Y_init = 20*sin(sqrt(X_init[,1]^2 + X_init[,2]^2))
  ran = matrix(c(0, 11), 2, 2)
  
  hpars = c(1.25, 1.25, 10, 
            1.1, 1.1, 0.1,
            0.1, 0.1)
  
  GP1 = CRNLHood$new(X_init, Y_init, ran, copula = F)
  GP1$Lhood_Prep()
  
  GP1$Refresh(Hpars = hpars, learnHpars=7, ymean=0)
  
  for(s in 0:2){
    x2  = matrix(c(4.25, 7.25, 1), 1)
    f   = function(x)GP1$kernel(matrix(c(x,s),1), x2)[1]
    df  = function(x)unlist(GP1$dkernel.dx1(matrix(c(x,s),1), x2))
    
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(2), 
                     h   = 0.001, 
                     tol = 0.01, 
                     msg = paste("CRN kernel x1, seed", s)
                     )
  }
}

test_CRN_post_kernel_grad = function(){
  
  # make some toy data
  set.seed(1)
  seeds = sample(1:10, replace = T, size = 100)
  
  X_init = cbind(rep(1:10, 10), rep(1:10, each=10), seeds)
  Y_init = 20*sin(sqrt(X_init[,1]^2 + X_init[,2]^2))
  ran = matrix(c(0, 11), 2, 2)
  
  hpars = c(1.25, 1.25, 10, 
            1.1, 1.1, 0.1,
            0.1, 0.1)
  
  GP1 = CRNLHood$new(X_init, Y_init, ran, copula = F)
  GP1$Lhood_Prep()
  
  GP1$Refresh(Hpars = hpars, learnHpars=7, ymean=0)
  
  x2  = matrix(c(4.25, 7.25, 1), 1)
  
  # posterior covraince with a fixed point
  for(s in 0:1){
    f   = function(x)as.numeric(GP1$COV(matrix(c(x,s),1), x2))
    df  = function(x)as.numeric(GP1$dCOV.dx1(matrix(c(x,s),1), x2))
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(2), 
                     h   = 0.001, 
                     tol = 0.05, 
                     msg = paste("CRN post kernel x1, seed", s)
    )
    
  }
  
  # posterior variance
  for(s in 0:1){
    f   = function(x)GP1$COV(matrix(c(x,s),1), matrix(c(x,s),1))[1]
    df  = function(x)unlist(GP1$dCOVxx.dx(matrix(c(x,s),1)))
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(2), 
                     h   = 0.001, 
                     tol = 0.05, 
                     msg = paste("CRN post kernel(x1, x1), seed", s)
    )
  }
  
  # posterior warping/shape changing
  for(s in 0:1){
    f   = function(x)GP1$SIGT_xr(matrix(c(x,s),1), x2)[1]
    df  = function(x)GP1$dSIGT.dx_xr(matrix(c(x,s),1), x2)[1,]
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(2),
                     h   = 0.00001, 
                     tol = 0.05, 
                     msg = paste("CRN post SIGT(x1; x2), seed", s)
    )
  }
  
  # posterior warping/shape changing
  for(s in 0:1){
    f   = function(x)GP1$SIGT_xr(matrix(c(x,s),1), matrix(c(x,0),1))[1]
    df  = function(x)GP1$dSIGT.dx_xr(matrix(c(x,s),1), matrix(c(x,s),1))[2,]
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(2),
                     h   = 0.00001, 
                     tol = 0.05, 
                     msg = paste("CRN post SIGT(x1; x1), seed", s)
    )
  }
  
  # posterior warping/shape changing
  f   = function(x)GP1$COV(matrix(c(x,1),1), matrix(c(x,0),1))[1]
  df  = function(x)unlist(GP1$dCOVxx0.dx(matrix(c(x,1),1)))
  tester_gradients(f, 
                   df, 
                   x0  = 10*runif(2),
                   h   = 0.00001, 
                   tol = 0.05, 
                   msg = "CRN post COV(x1, x1_0)"
  )

}


##############################
# acquisition function tests

test_KGCB_expectation = function(){
  
  cat("\nTesting KGCB expectation")
  
  set.seed(1)
  
  KGv=sapply(1:10, function(i){
    a = runif(5)
    b = 0.5 - runif(5)
    emp = mean(apply(a + outer(b, rnorm(100000)), 2, max)) - max(a)
    thr = KGCBfilter_grad(a, b)$KG
    c(emp, thr)
  })
  
  cat("\nMC :", KGv[1,])
  cat("\nKG :", KGv[2,], "\n")
  assert("MC and thry diff is within 10%", max(abs((KGv[1,]-KGv[2,])/KGv[1,]))<0.3)
  cat("\nTest Passed!\n")
  
}

test_KGCB_grad = function(){
  
  
  KGCB = function(ab){ dd = length(ab)/2; KGCBfilter_grad(ab[1:dd], ab[-(1:dd)])$KG }
  
  dKGCB = function(ab){ dd = length(ab)/2; dKG=KGCBfilter_grad(ab[1:dd], ab[-(1:dd)]); c(dKG$dmu, dKG$dsig) }
  
  set.seed(1)
  a = runif(5)
  b = 0.5 - runif(5)
  
  tester_gradients(f   = KGCB,
                   df  = dKGCB,
                   x0  = c(a,b),
                   h   = 0.0001,
                   msg = "KGCB(a,b) test 1")
  
  set.seed(2)
  a = runif(5)
  b = 0.5 - runif(5)
  cat("\n")
  cat("\n")
  tester_gradients(f   = KGCB,
                   df  = dKGCB,
                   x0  = c(a,b),
                   h   = 0.0001,
                   msg = "KGCB(a,b) test 2")
}

test_KG_gradx = function(){
  
  
  # make some toy data
  set.seed(1)
  X1    = sample(1:10, replace = T, size = 10)
  X2    = sample(1:10, replace = T, size = 10)
  seeds = sample(1:4,  replace = T, size = 10)
  
  X_init = cbind(X1, X2, seeds)
  Y_init = 20*sin(sqrt(X_init[,1]^2 + X_init[,2]^2))
  ran = matrix(c(0, 11), 2, 2)
  
  hpars = c(2, 2, 10, 
            2, 2, 0.1,
            0.1, 0.1)
  
  GP1 = CRNLHood$new(X_init, Y_init, ran, copula = F)
  GP1$Lhood_Prep()
  
  GP1$Refresh(Hpars = hpars, learnHpars=7, ymean=0)
  
  # make the KG function
  Xr = cbind(LHCran(117, ran), 0)
  plot(X_init, cex=(Y_init-min(Y_init)+40)/20, pch=19)
  points(Xr)
  
  CRNKG = Make_CRNKG_grad(GP1, Xr)
  
  for(i in 1:11){
    f = function(x)CRNKG$KG(c(x,i))
    df = function(x)CRNKG$dKG(matrix(c(x, i), 1))$dKG
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(2), 
                     h   = 0.001, 
                     tol = 0.1, 
                     msg = paste("dCRKG/dx seed", i)
                     )
  }
  
}

test_PWKG_gradx = function(){
  
  
  # make some toy data
  set.seed(1)
  X1    = sample(1:10, replace = T, size = 10)
  X2    = sample(1:10, replace = T, size = 10)
  seeds = sample(1:4,  replace = T, size = 10)
  
  X_init = cbind(X1, X2, seeds)
  Y_init = 20*sin(sqrt(X_init[,1]^2 + X_init[,2]^2))
  ran = matrix(c(0, 11), 2, 2)
  
  hpars = c(2, 2, 10, 
            2, 2, 0.1,
            0.1, 0.1)
  
  GP1 = CRNLHood$new(X_init, Y_init, ran, copula = F)
  GP1$Lhood_Prep()
  
  GP1$Refresh(Hpars = hpars, learnHpars=7, ymean=0)
  
  # make the KG function
  Xr = cbind(LHCran(117, ran), 0)
  plot(X_init, cex=(Y_init-min(Y_init)+40)/20, pch=19)
  points(Xr)
  
  PWKG = Make_PWKG_grad(GP1, Xr)
  
  for(i in 1:11){
    f = PWKG$PWKG
    df = PWKG$dPWKG
    tester_gradients(f, 
                     df, 
                     x0  = 10*runif(4), 
                     h   = 0.001, 
                     tol = 0.1, 
                     msg = paste("dPWKG/dx point", i)
    )
  }
  
}

##############################
# Test gradient ascender

test_optimizer_KG = function(){
  
  
  # make some toy data
  set.seed(1)
  X1    = sample(1:10, replace = T, size = 10)
  X2    = sample(1:10, replace = T, size = 10)
  seeds = sample(1:4,  replace = T, size = 10)
  
  X_init = cbind(X1, X2, seeds)
  Y_init = 20*sin(sqrt(X_init[,1]^2 + X_init[,2]^2))
  ran = matrix(c(0, 11), 2, 2)
  
  hpars = c(2, 2, 10, 
            2, 2, 0.1,
            0.1, 0.1)
  
  GP1 = CRNLHood$new(X_init, Y_init, ran, copula = F)
  GP1$Lhood_Prep()
  
  GP1$Refresh(Hpars = hpars, learnHpars=7, ymean=0)
  
  # make the KG function
  Xr = cbind(LHCran(117, ran), 0)
  plot(X_init, cex=(Y_init-min(Y_init)+40)/20, pch=19)
  points(Xr)
  
  CRNKG = Make_CRNKG_grad(GP1, Xr)
  
  f = function(x)CRNKG$KG(c(x,11))
  df = function(x)CRNKG$dKG(matrix(c(x, 11), 1))$dKG
  
  for (i in 1:10){
    set.seed(i)
    output1 = Optimizer_2(f, df, ran, N0 = 10, Na = 0)

    set.seed(i)
    output2 = Optimizer_2(f, df, ran, N0 = 10, Na = 6, maxevals = 100)

    cat("\n", output1$fmax, output2$fmax, "\n")
    assert("grad ascent must work", output2$fmax>=output1$fmax)
  }
  
  for (i in 1:10){
    set.seed(i)
    x0 = 10*runif(2)
    f0 = f(x0)
    
    output1 = Optimizer_2(f, df, ran, N0 = 0, Na = 0, x0=x0, maxevals = 200)
    
    cat("\n", f0, output1$fmax, "\n")
    assert("single initial point", output1$fmax>=f0)
  }
  
  cat("\nTest Passed!\n")
}

# test_KG_optimizer()


if(T){
  AA = ls()
  tests = AA[grepl("test_", AA)]
  sapply(tests, function(t)match.fun(t)())
}



