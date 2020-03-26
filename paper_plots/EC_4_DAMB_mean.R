# setwd("/Users/pearce/CRN/")
source("plotting/plot_utils.R")

Results = readRDS("/home/michael/CopiedResults/8001")


M = sapply(Results, function(r)r$method)
R = Get_multiple_OC(Results, M, N0=21, N1=1000)

S = Get_multiple_seed_reuse(Results, M, N0=23, N1=1000)

Tim = Get_multiple_Timings(Results, M, N0=21, N1=1000)

HP_7 = Get_multiple_Hypers(Results, M, 7)
HP_14 = Get_multiple_Hypers(Results, M, 14, N0=21, N1=1000)
HP_15 = Get_multiple_Hypers(Results, M, 15, N0=21, N1=1000)
HP_16 = Get_multiple_Hypers(Results, M, 16, N0=21, N1=1000)

cols = rep("black",10)

"
Experiment 8001
- ambulance, 1800s per sim, all results, KG, PW, CRN, 400 runs each, upto 1000 samples
- objective: mean journey time
- discretization is old obs and LHC
- CRN KG optimizer evaluates all seeds at the final x.
- PWKG optimizer does KG first, then does pair with first frozen then pair with both free
- method numbers are a little fucked up in the results
"

# M values
# 2 KG
# 5 PWKG CS
# 6 broken BS
# 7 CRNKG CSW
# 8 CRNKG CS
# 9 CRNKG CSW

R = lapply(R, function(r){r[,2] = abs(r[,2]); r})

cols[1] = "blue" # KG
cols[2] = "red" # PWKG CS
cols[3] = "green" # faulty
cols[4] = "darkred" # PWKG+W
cols[5] = "darkgreen" # CRNKG PW
cols[6] = "green" #CRNKG CSW


pdf("~/Dropbox/PhD/CRN/CRNKG paper/Operations-Research-Paper/April2019/pics/Code_for_plots/EC_4_DAMB_mean.pdf",
    height = 5, width = 8)
par(mfrow=c(2,3),mar = c(4.1,3.1,3.1,1.1))

###################################################################
###################################################################
######################### TOP ROW #################################

# FIRST PLOT PROFIT
Res_plot(R, cols = cols,  erscale = 0.5, yr=c(0.142, 0.155))
title(main = "Average Journey Time", font.main=1, xlab="N", ylab="")
legend(legend = c("KG", "KG_PW", "KG-PW-Bias", "KG-CRN-CS", "KG-CRN"),
       col = c("blue", "green", "darkgreen", "red", "darkred"),
       lty=rep(1,5),
       lwd=4,
       x="bottomleft")

# SECOND PLOT
Res_plot(S, cols = cols)
lines(c(20, 1000), c(0.5, 0.5), lty=2)
title(main="Seed Reuse Frequency", xlab="N", font.main=1)

# THIRD PLOT CUM TIME
Res_plot(Tim, cols = cols)
# lines(c(20, 500), c(0.5, 0.5), lty=2)
title(main="Cumulative Time, hours", xlab="N", font.main=1)


###################################################################
###################################################################
#######################  BOTTOM ROW  ##############################

# first PLOT BIAS PARAMETER
Res_plot(HP_14, cols=cols, yr=c(0, 7e-4))
title(main=expression(paste("Bias Parameter,", {sigma[b]}^2)),
      xlab="N",
      ylab=expression(paste("Profit"^2)))

# THIRD PLOT OFFSET PARAMETER
Res_plot(HP_15, cols=cols, yr=c(0, 7e-4))
title(main=expression(paste("Offset Parameter,", {eta}^2)),
      xlab="N",
      ylab=expression(paste("Profit"^2)))

# FOURTH PLOT WHITE NOISE PARAMETER
Res_plot(HP_16, cols=cols, yr=c(0, 7e-4))
title(main=expression(paste("White Noise Parameter,", {sigma[w]}^2)),
      xlab="N",
      ylab=expression(paste("Profit"^2)))


dev.off()



