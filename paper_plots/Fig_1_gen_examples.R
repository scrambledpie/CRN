# .rs.restartR()
# Sys.sleep(5)
rm(list = ls())
setwd("~/Dropbox/PhD/CRN/CRN/")
source("CRN_BO/Acquisition_funs.R")


Hypers = c(10, 1^2, 10, 0.1, 0.3, 0.05)
################################################################################
# Make test data
Make_testFun = function(TruePars = Hypers, seed=9){
  set.seed(seed)
  # TruePars = c(5, 100^2, 5, wig, off, wit)
  Nseeds   = 115
  X_domain = (1:1000)*100/1000
  
  
  SEkernel = function(x1,x2) TruePars[2]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[1]^2)
  X_f = X_domain
  NX = length(X_domain)
  Y_f = mvrnorm(1, rep(0, NX), SEkernel(X_f, X_f))
  Y_f = matrix(Y_f, NX, Nseeds+1)
  Y_f = Y_f - mean(Y_f)
  
  SEkernelw = function(x1,x2)TruePars[4]*exp(-0.5*outer(x1,x2,'-')^2/TruePars[3]^2)
  Y_wiggles = t(mvrnorm(Nseeds, rep(0, NX), SEkernelw(X_f,X_f)))
  Y_wiggles = cbind(0, Y_wiggles)
  
  Y_Offsets = c(0, rnorm(Nseeds, 0, sqrt(TruePars[5])))
  Y_Offsets = matrix(Y_Offsets , NX, Nseeds+1, byrow = T)
  
  Y_white   = matrix(rnorm(NX*Nseeds, 0, sqrt(TruePars[6])), NX, Nseeds)
  Y_white   = cbind(0, Y_white)
  
  TestFunY  = Y_f + Y_wiggles + Y_Offsets + Y_white
  
  TestFunY2  = Y_f + Y_wiggles + Y_Offsets

  
  TestFun_i = function(xs){
    x_i = which.min(abs(xs[1]-X_domain))
    TestFunY[x_i, xs[2]+1]
  }
  
  TestFunNW_i = function(xs){
    x_i = which.min(abs(xs[1]-X_domain))
    TestFunY2[x_i, xs[2]+1]
  }
  
  TestFun = function(xs, nowhite=F){
    xs = matrix(xs, ncol=2)
    if(nowhite){
      apply(xs,1,TestFunNW_i)
    }else{
      apply(xs,1,TestFun_i)
    }
  }
  
  return(TestFun)
}

N0 = 4
N1 = 8
TestFun = Make_testFun()
xs = cbind(round((1:N0-0.5)*100/N0), rep(c(1,2), len=N0))
y = TestFun(xs)


par(mfcol=c(2,2), mar = c(4.1,5.1,3.1,1.1))


plot_Gen = function(X_r,
                    TruePars=Hypers, 
                    seed=2, 
                    truth_lines=F,
                    truth_points=F,
                    seed_lines=F,
                    seed_points=F,
                    Svals = 0:2,
                    yran=c(-2.3,3.7),
                    main="",
                    dense=F,
                    ...
                    ){
  
  
  TestFun = Make_testFun(TruePars, seed)
  
  # X_r = round(X_r)
  
  Y = sapply(Svals, function(s)TestFun(cbind(X_r, s)))
  
  Ynw = sapply(Svals, function(s)TestFun(cbind(X_r, s), nowhite=T))

  plot(c(0, 100), yran, col="white", xlab="X", ylab="Y", main=main, font.main=1)
  
  
  
  
  litcols = c("black", "salmon", "lightgreen")
  dkcols  =c("black", "red", "green")
  if(seed_points)for(s in 2:ncol(Y)){
    X_r1 = c(X_r, X_r)
    if(dense){Y2 = Ynw[,s]}else{Y2=Y[,s]}
    
    Y_s1 = c(Y[,s], Y2)
    # points(X_r1, Y_s1, col=1, pch=19, cex=1.35)
    # if(dense)lines(X_r, Ynw[,s], col="black", lwd=13)
    points(X_r1, Y_s1, col="white", pch=19, cex=1)
    points(X_r1, Y_s1, col=litcols[s], pch=19, cex=1)
    if(dense)lines(X_r, Ynw[,s], col=litcols[s], lwd=11)
    
    
  }
  if(seed_lines)for(s in 2:ncol(Y)){
    # lines(X_r, Ynw[,s], col=1, lwd=4)
    lines(X_r, Ynw[,s], col=dkcols[s], lwd=2)
    # topY = which.max(Ynw[,s])
    # points(X_r[topY], Ynw[topY,s], pch=8, cex=2, col=s)  
  }
  
  if(truth_points)points(X_r, Y[,1], col=1, pch=19)
  if(truth_lines) lines(X_r, Y[,1], col=1, lwd=4)
  
  

}


pdf("~/Dropbox/PhD/CRN/CRNKG paper/Operations-Research-Paper/April2019/pics/Code_for_plots/Fig_1_gen_examples.pdf",
    height = 6, width = 8)
par(mfcol=c(2,2), mar = c(4.1,4.1,3.1,1.1))

seed=71
# Just a few points
N0 = 25
X_r = (1:N0-0.5)*100/N0
plot_Gen(X_r, TruePars = Hypers, seed=seed, T,F,T,T, main="General Model", main.font=1)
legend(10, 3.7, lwd=c(2,2,2), pch=c(NA, 19,19), col=c("black", "red", "green"), legend=c("Target", "s=1", "s=2"), ncol=3)

# Wiggle only
N0 = 25
X_r = (1:N0-0.5)*100/N0
Hypers_wig = Hypers
Hypers_wig[4] = sum(Hypers[4:6])
Hypers_wig[5] = 0
Hypers_wig[6] = 0
plot_Gen(X_r, TruePars = Hypers_wig, seed=seed, T,F,T,T, main="Biases Only", yran=c(-3,3))


# full coreraltion
N0 = 25
X_r = (1:N0-0.5)*100/N0
Hypers_full = Hypers
Hypers_full[4] = 0
Hypers_full[5] = sum(Hypers[4:6])
Hypers_full[6] = 0
plot_Gen(X_r, TruePars = Hypers_full, seed=seed, T, F, T, T, main="Offsets Only")


# Dense X
N0=500
X_r = 1:N0*100/N0 
Hypers_dense = Hypers
Hypers_dense[4] = 0
Hypers_dense[5] = sum(Hypers[4:5])
plot_Gen(X_r, TruePars = Hypers_dense, seed=seed, T, F, T, T, main="Offsets and White Noise, Dense X", yran=c(-3,3), dense = T)
dev.off()




