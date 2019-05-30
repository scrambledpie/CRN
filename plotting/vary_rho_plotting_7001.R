source("plotting/plot_utils.R")

Results = readRDS("~/CopiedResults/7001")

# Get the details of each element in the Results array
Lens = sapply(Results, function(r)max(r$Cost$N))
# RHOS = sapply(Results, function(r)sum(tail(r$Hpars,1)[[1]][4:6]))
RHOS = sapply(Results, function(r)sum(tail(r$Hpars,1)[[1]][4]))
M    = sapply(Results, function(r)r$method)
uRHOS = sort(unique(RHOS))
nR = length(uRHOS)

# now make a vector that is unique for each setting
MR = sapply(1:length(RHOS), function(i){ 
    ri = which(RHOS[i]==uRHOS)
    mi = which(M[i]==c(2,4,5))
    Xid = ri + nR*(mi-1)
})

# new vectors with properties of each experiment
MR_cols = rep(c("blue", "red", "green"), each=nR)
MR_MR   = sort(unique(MR))
# MR_M    = rep(c(2, 4, 5), each=nR)
# MR_R    = rep(uRHOS, 3)
MR_M    = c(2, rep(c(4,5), each=nR))
MR_R    = c(2, rep(uRHOS, 2))
MR_OC   = Get_multiple_OC(Results, MR, with_max = T)
MR_SF   = Get_multiple_seed_reuse(Results, MR)

# plot all the OC curves, one for each rho
par(mfcol=c(2,10))
ran = Get_xy_ran(MR_OC)
for( rho in uRHOS[-1]){
  SS = c(T, rho==MR_R[-1])
  # SS = rho==MR_R
  
  Res_plot(MR_OC[SS], ran$xr, c(1,ran$yr[2]), cols=c("blue", "red", "green"), main=signif(rho/2500,1), log="y")
  Res_plot(MR_SF[SS], ran$xr, c(0, 1), cols=c("blue", "red", "green"), main=signif(rho/2500,1))
  lines(ran$xr, c(0.5,0.5), lty=2)
  
  
}


NN = 96
par(mfrow=c(1,1))
plot_errors(MR_OC, MR_M, MR_R, NN, cols=c("blue", "red", "green"))

KG_f = MR_OC[[1]][NN,][2:3]
KG_f = cbind(uRHOS, matrix(KG_f, nR, 2, byrow = T))
Plot_pnts_bars(KG_f, "blue")
title(main="100 samples, KG=blue, CRN=red, PW=Gr", xlab = "wiggle ratio", ylab = "OC")



