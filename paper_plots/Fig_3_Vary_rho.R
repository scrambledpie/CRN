source("plotting/plot_utils.R")

Results = readRDS("~/CopiedResults/1067")

# Res# Get the details of each element in the Results array
Lens = sapply(Results, function(r)max(r$Cost$N))
RHOS = sapply(Results, function(r)r$rho)
M    = sapply(Results, function(r)r$method)

MR = RHOS*10 + M*100

uMR = sort(unique(MR))
uRHOS = sort(unique(RHOS))

# new vectors with properties of each experiment

MR_OC   = Get_multiple_OC(Results, MR, with_max = T)
MR_SF   = Get_multiple_seed_reuse(Results, MR)
MR_NS   = Get_multiple_novel_seed_count(Results, MR)


strip = function(res){
  res[1,2]=0
  res[47,2] = res[46,2]
  res
}

pdf("~/Dropbox/PhD/CRN/CRNKG paper/Operations-Research-Paper/April2019/pics/Code_for_plots/Fig_3_Vary_rho.pdf",
    height = 6,
    width = 8)
par(mfrow=c(2,2), mar = c(4.1,4.1,3.1,1.1))


# FIRST PLOT OC for best rho
SS = abs((uMR%%100)-10)<0.0001
Res_plot(MR_OC[SS], ran$xr, c(0.00001,30), cols=c("blue", "red", "green"), log="")
title(xlab="N", ylab="OC", 
      main="Opportunity Cost for Full Correlation",
      font.main=1)
legend(
  col    = c("blue", "green", "red"),
  legend = c("KG", "KG-PW", "KG-CRN"),
  lwd = rep(2,3),
  x = 20, y=30
)
# lines(c(25,25), c(0, 25), lty=2)

# SECOND PLOT OC  VS RHO
plot_errors(MR_OC, floor(uMR/100)/10, rep(0:10, 3)/10, 45, c("blue", "red", "green"))
title(main="OC after 50", xlab=expression(rho), ylab="OC", font.main=1)


# THIRD PLOT SEED REUSE HIGH+LOW
MR_SF[uMR>499] = lapply(MR_SF[uMR>499], strip)
Res_plot(MR_SF[SS], 
         c(5, 50), 
         c(0, 1), 
         cols=c("blue", "red", "green"), 
         ylab="Frequency")
title(main="Seed Reuse Frequency", font.main=1) 
legend(
  col    = c("black", "black"),
  legend = c(expression(paste(rho, "= 1")), expression(paste(rho, "= 0"))),
  lwd = rep(2,3),
  lty = c(1, 2),
  x = 5, y=0.8
)

Plot_mean_er(MR_SF[uMR==500][[1]], "green", lty=2)
lines(MR_SF[uMR==500][[1]][,1], MR_SF[uMR==500][[1]][,2], col="darkgreen", lty=2)
lines(MR_SF[uMR==510][[1]][,1], MR_SF[uMR==510][[1]][,2], col="darkgreen", lty=1)
# lines(MR_SF[uMR==500][[1]][,1], MR_SF[uMR==500][[1]][,2], col="darkgreen", lty=2)
Plot_mean_er(MR_SF[uMR==400][[1]], "red", lty=2)
lines(MR_SF[uMR==400][[1]][,1], MR_SF[uMR==400][[1]][,2], col="darkred", lty=2)
lines(MR_SF[uMR==410][[1]][,1], MR_SF[uMR==410][[1]][,2], col="darkred", lty=1)
lines(c(5,50), c(0.5,0.5), lty=3)
title(xlab="N")

# FOURTH PLOT 
plot_errors(MR_SF, floor(uMR/100)/10, rep(0:10, 3)/10, 45, c("blue", "red", "green"))
title(main="Seed Reuse after 50", xlab=expression(rho), ylab="Frequency", font.main=1)
lines(c(0,1), c(0.5,0.5), lty=3)


dev.off()


