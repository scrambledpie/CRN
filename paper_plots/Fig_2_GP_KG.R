# .rs.restartR()
# Sys.sleep(5)
rm(list = ls())
setwd("~/Dropbox/PhD/CRN/CRN/")
source("CRN_BO/Acquisition_funs.R")


Hypers = c(10, 1^2, 10, 0.2, 0.15, 0.2)
################################################################################
# Make test data
Make_testFun = function(TruePars = Hypers, seed=8){
  set.seed(seed)
  # TruePars = c(5, 100^2, 5, wig, off, wit)
  Nseeds   = 115
  X_domain = 1:100
  
  SEkernel = function(x1,x2) TruePars[2]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[1]^2)
  X_f = X_domain
  Y_f = mvrnorm(1, rep(0,length(X_f)), SEkernel(X_f, X_f))
  Y_f = matrix(Y_f, 100, Nseeds+1)
  Y_f = Y_f - mean(Y_f)
  
  SEkernelw = function(x1,x2)TruePars[4]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[3]^2)
  Y_wiggles = t(mvrnorm(Nseeds, rep(0, length(X_f)), SEkernelw(X_f,X_f)))
  Y_wiggles = cbind(0, Y_wiggles)
  
  Y_Offsets = c(0, rnorm(Nseeds, 0, sqrt(TruePars[5])))
  Y_Offsets = matrix(Y_Offsets , 100, Nseeds+1, byrow = T)
  
  Y_white   = matrix(rnorm(100*Nseeds, 0, sqrt(TruePars[6])), 100, Nseeds)
  Y_white   = cbind(0, Y_white)
  
  TestFunY  = Y_f + Y_wiggles + Y_Offsets + Y_white
  
  TestFunY2  = Y_f + Y_wiggles + Y_Offsets
  
  TestFun = function(xs, nowhite=F){
    if (any(xs != round(xs))) stop("give integer inputs")
    xs = matrix(xs, ncol=2)
    if(nowhite){
      apply(xs,1,function(xs) TestFunY2[xs[1], xs[2]+1])
    }else{
      apply(xs,1,function(xs) TestFunY[xs[1], xs[2]+1])
    }
  }
  
  return(TestFun)
}

N0 = 4
N1 = 8
TestFun = Make_testFun()
xs = cbind(round((1:N0-0.5)*100/N0), rep(c(1,2), len=N0))
y = TestFun(xs)


plot_GP_and_KG = function(xs, 
                          y, 
                          TruePars = Hypers, 
                          plotting=T, 
                          KGmax=0.22, 
                          legends=T, 
                          legendx = 0, 
                          gpran=c(-1.8,2.0)){
  
  
  stopifnot(
    all(xs>=1),
    all(xs<=100)
  )
  
  GP1 = CRNLHood$new(xs,
                     y,
                     matrix(c(0,100), 2,1), 
                     copula=F
  )
  
  GP1$Refresh(Hpars = TruePars, learnHpars = 7)
  
  
  X_r = cbind(1:100, 0)
  
  KG1 = Make_CRNKG_grad(GP1, X_r)$KG
  
  KG = function(xs){xs=matrix(xs, ncol=2); apply(xs, 1, KG1)}
  
  if(plotting){
    # FIRST PLOT GP MODEL
    plot(c(0,100), 
         gpran,
         col="white",
         xlab="X",
         ylab="Y",
         main=paste("GP Model, n=", nrow(xs), sep=''),
         font.main=1
    )
    
    if(legends)legend(x="topright", 
                      legend=c("truth","s=0", "s=1", "s=2"),
                      col=c("grey",1:3),
                      lty=c(2,1,1,1),
                      lwd=2,
                      ncol=1
    )
    
    Y_true = TestFun(cbind(X_r[,1], 0), nowhite=T)
    lines(X_r[,1], Y_true, col='grey', lty=2)
    
    for( s in 0:max(xs[,2])){
      xs_s = cbind(xs[xs[,2]==s,1]-0.1, s)
      Y_s0 = GP1$MU(xs_s)[,1]
      Y_s = y[xs[,2]==s]
      
      Y_r = GP1$MU(cbind(X_r[,1], s))[,1]
      
      plot_Y = c(Y_r, Y_s0, Y_s, Y_s0)
      plot_X = c(X_r[,1], xs_s[,1]-0.001, xs_s[,1], xs_s[,1]+0.001)
      
      plot_order = order(plot_X)
      plot_X = plot_X[plot_order]
      plot_Y = plot_Y[plot_order]
      
      lines(plot_X, plot_Y, col=s+1)
      if(s==0)lines(X_r[,1], Y_r, col=s+1, lwd=2)
      
      if(sum(xs[,2]==s)>0){
        keep = xs[,2]==s
        points(xs[keep, 1], y[keep], col=s+1, pch=19)
      }
    }
    
    
    # SECOND PLOT KG VALUES
    if (nrow(xs)==4)main_tit=expression(paste("KG"^CRN,", n=4", sep=""))
    if (nrow(xs)==5)main_tit=expression(paste("KG"^CRN,", n=5", sep=""))
    if (nrow(xs)==6)main_tit=expression(paste("KG"^CRN,", n=6", sep=""))
    if (nrow(xs)==7)main_tit=expression(paste("KG"^CRN,", n=7", sep=""))
    if (nrow(xs)==8)main_tit=expression(paste("KG"^CRN,", n=8", sep=""))
    if (nrow(xs)==9)main_tit=expression(paste("KG"^CRN,", n=9", sep=""))
    if (nrow(xs)==10)main_tit=expression(paste("KG"^CRN,", n=10", sep=""))
    plot(c(0,100),
         c(0, KGmax),
         col="white",
         xlab="X",
         ylab=expression(paste("KG"^CRN,"(x,s)", sep="")),
         main=main_tit
    )
    
    if(legends)legend(x="topright",
                      legend=c("max  ", "s=1", "s=2", "s=3"),
                      col=1:4,
                      lwd=c(NA,2,2,2),
                      pch=c(8,NA,NA,NA),
                      ncol=1)
  }
  
  topKG = -Inf
  topxs = 0
  for( s in 1:max(xs[,2]+1) ){
    X_r_s = sort(c(X_r[,1], xs[,1]-0.1, xs[,1]+0.1))
    Y_s = KG(cbind(X_r_s, s))
    if(plotting) lines(X_r_s, Y_s, col=s+1)
    
    if (max(Y_s)>topKG){
      best_i = which.max(Y_s)
      topKG = Y_s[best_i]
      topxs = c(X_r_s[best_i], s)
    }
    
    
    # if(sum(xs[,2]==s)>0&plotting){
    #   keep = xs[,2]==s
    #   points(xs[keep, 1], rep(0, sum(keep)), col=s+1, pch=19)
    # }
    
  }
  
  if(plotting) points(topxs[1], topKG, pch=8)
  
  return(topxs)
}


pdf("~/Dropbox/PhD/CRN/CRNKG paper/Operations-Research-Paper/April2019/pics/Code_for_plots/Fig_2_GP_KG.pdf",
    height = 6,
    width = 8)
par(mfcol=c(2,2), mar = c(4.1,4.3,3.1,1.1))

topxs = plot_GP_and_KG(xs, y, plotting = T)

while(nrow(xs)<N1){
  topxs = plot_GP_and_KG(xs, y, plotting = F)
  xs = rbind(xs, topxs)
  print(topxs)
  y = c(y, TestFun(topxs))
}

topxs = plot_GP_and_KG(xs, y, plotting = T, legends = F)

dev.off()