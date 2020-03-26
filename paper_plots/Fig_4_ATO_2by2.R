# setwd("/Users/pearce/CRN/")
setwd("~/Dropbox/PhD/CRN/CRN/")
source("plotting/plot_utils.R")

Results = readRDS("/home/michael/CopiedResults/5006")


M = sapply(Results, function(r)r$method)
R = Get_multiple_OC(Results, M, N0=21, N1=500)

S = Get_multiple_seed_reuse(Results, M, N0=23, N1=500)

HP_7 = Get_multiple_Hypers(Results, M, 9)
HP_14 = Get_multiple_Hypers(Results, M, 18, N0=21, N1=500)
HP_15 = Get_multiple_Hypers(Results, M, 19, N0=21, N1=500)
HP_16 = Get_multiple_Hypers(Results, M, 20, N0=21, N1=500)

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

pdf("~/Dropbox/PhD/CRN/CRNKG paper/Operations-Research-Paper/April2019/pics/Code_for_plots/Fig_4_ATO_2by2.pdf",
    height = 6, width = 8)
par(mfrow=c(2,2), mar = c(4.1,4.1,3.1,1.1))

# FIRST ROW
Res_plot(R, cols = cols, erscale = 0.5, xlab="N")
title(main = "Average Profit", font.main=1, ylab="Profit")
# Res_plot(HP_14, cols=cols, main="Bias Function", yr=c(0, 1200))
legend(legend = c("KG", "KG-PW", "KG-PW-bias", "KG-CRN-CS", "KG-CRN"),
       col = c("blue", "green", "darkgreen", "red", "darkred"),
       lty=rep(1,5),
       lwd=4,
       x="bottomright")
# Res_plot(HP_15, cols=cols, main="Offset", yr=c(0, 1200))
# Res_plot(HP_16, cols=cols, main="White Noise", yr=c(0, 1200))
# 

# SECOND SEED REUSE
Res_plot(S, cols = cols, xlab="N")
title(main = "Seed Reuse", font.main=1, ylab="Frequency")
lines(c(20, 500), c(0.5, 0.5), lty=2)

CRNbars = Get_seed_bars(Results[M==4], 6)
barplot(CRNbars, names.arg = 1:ncol(CRNbars), ylim=c(0, 500),
        xlab="Seed",
        ylab="Samples")
title(main=expression(paste("KG"^CRN, "-CS Seed Allocation", sep="")))

CRNbars = Get_seed_bars(Results[M==6], 6)
barplot(CRNbars, names.arg = 1:ncol(CRNbars), ylim=c(0, 500),
        xlab="Seed",
        ylab="Samples")
title(main=expression(paste("KG"^CRN, " Seed Allocation", sep="")))
dev.off()
# 
# CS_PWbars = Get_seed_pairs_singles(Results[M==5])
# CSW_PWbars = Get_seed_pairs_singles(Results[M==7])
# barplot(c(CS_PWbars, CSW_PWbars), main="PW-KG Seeds", 
#         names.arg=c("CS\nsingles", "pairs", "Bias\nsingles", "pairs"),
#         las=2)



# Save as PDF 8x4 inches


if(F){
  
  means = sapply(R, function(r)r[480,2])
  ers = (sapply(R, function(r)r[480,3])/2)^2
  best = which.max(means)
  pnorm( (means[best] - means)/sqrt(ers[best]+ers) ) <0.95
  
}
