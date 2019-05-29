source("plotting/plot_utils.R")

Results = readRDS("~/CopiedResults/7001")

# Get the details of each element in the Results array
Lens = sapply(Results, function(r)max(r$Cost$N))
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
MR_M    = rep(c(2,4,5), each=nR)
MR_R    = rep(uRHOS, 3)
MR_OC   = Get_multiple_OC(Results, MR, with_max = T)
MR_SF   = Get_multiple_seed_reuse(Results, MR)

# plot all the OC curves, one for each rho
par(mfrow=c(2,1))
ran = Get_xy_ran(MR_OC)
for( rho in uRHOS){
  Res_plot(MR_OC[rho==MR_R], ran$xr, c(1e-2,ran$yr[2]), cols=c("blue", "red", "green"), main=rho, log="y")
  
  Res_plot(MR_SF[rho==MR_R], ran$xr, c(0,1), cols=c("blue", "red", "green"), main=rho)
  lines(ran$xr, c(0.5,0.5), lty=2)
  
  
}






