# setwd("/Users/pearce/CRN/")

# source("~/Dropbox/PhD/CRN/CRN/plotting/plot_utils.R")
# Results = readRDS("/home/michael/CopiedResults/8001")

# BO algorithm results
source("/Volumes/DataStorage/Dropbox/PhD/CRN/CRN/plotting/plot_utils.R")
Results = readRDS("/Volumes/DataStorage/CopiedResults/8001")

# Add on the SOSA result
Results2 = readRDS("/Volumes/DataStorage/CopiedResults/113_CRN")
Results2 = lapply(Results2, function(r){r$method=10; r})

Results3 = c(Results, Results2)
M3 = sapply(Results3, function(r)r$method)

R = Get_multiple_OC(Results3, M3, N0=21, N1=1000)

M = sapply(Results, function(r)r$method)

S = Get_multiple_seed_reuse(Results, M, N0=23, N1=1000)

HP_7 = Get_multiple_Hypers(Results, M, 7)
HP_14 = Get_multiple_Hypers(Results, M, 14, N0=21, N1=1000)
HP_16 = Get_multiple_Hypers(Results, M, 16, N0=21, N1=1000)
HP_15 = Get_multiple_Hypers(Results, M, 15, N0=21, N1=1000, med=TRUE)

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
cols[6] = "purple" #CRNKG CSW


cols = c(cols, "purple", NULL, NULL, "purple")

savepdf=F
if(savepdf)pdf("~/Dropbox/PhD/CRN/CRNKG paper/Operations-Research-Paper/April2019/pics/Code_for_plots/Fig_5_DAMB_2by2.pdf",
               height = 6, width = 8)

par(mfrow=c(2,2), mar = c(4.1,4.1,3.1,1.1))

# FIRST ROW
# Res_plot(R, cols=cols, main = "", yr=c(0.142, 0.154), erscale = 0.5,xlab = "N", ylab = "Time")
Res_plot(R, cols=cols, main = "", yr=c(0.142, 0.167), erscale = 0.5,xlab = "N", ylab = "Time")
title(main="Average Time to Patients", font.main=1)
# Res_plot(HP_14, cols=cols, main="Bias Function", yr = c(0, 0.0006))
leg_names = c(expression(paste('SOSA')),
              expression(paste('KG')),
              expression(paste('KG'^PW)),
              expression(paste('KG'^PW,'-bias')),
              expression(paste('KG'^CRN,'-CS')),
              expression(paste('KG'^CRN)))
leg_names = c(expression(paste('SOSA')),
              expression(paste('KG')),
              expression(paste('KG-PW')),
              expression(paste('KG-PW-bias')),
              expression(paste('KG-CRN-CS')),
              expression(paste('KG-CRN')))
legend(legend = leg_names,
       col = c("purple", "blue", "green", "darkgreen", "red", "darkred"),
       lty=rep(1,5),
       lwd=rep(4, 5),
       x="topright")
# Res_plot(HP_15, cols=cols, main="Offset", yr = c(0, 0.0006))
# Res_plot(HP_16, cols=cols, main="White Noise", yr = c(0, 0.0006))


Res_plot(S, cols = cols, main = "", xlab="N")
lines(c(20, 1000), c(0.5, 0.5), lty=2)
title(main="Seed Reuse", font.main=1)

# SECOND ROW
CRNbars = Get_seed_bars(Results[M==4], 10)
barplot(CRNbars, names.arg = 1:ncol(CRNbars),
        main="",
        ylim=c(0, 600),
        xlab = "Seed", ylab="Samples")
title(main=expression(paste('KG'^CRN,'-CS Seed Allocation')), font.main=1)

CRNbars = Get_seed_bars(Results[M==6], 30)
barplot(CRNbars, names.arg = 1:ncol(CRNbars),
        main="",
        ylim=c(0, 300),
        xlab = "Seed", ylab="Samples")
title(main=expression(paste('KG'^CRN,' Seed Allocation')), font.main=1)

if(savepdf)dev.off()
# 
# CS_PWbars = Get_seed_pairs_singles(Results[M==5], topseed = 1000)
# CSW_PWbars = Get_seed_pairs_singles(Results[M==7], topseed = 1000)
# barplot(c(CS_PWbars, CSW_PWbars), main="PW-KG Seeds", 
#         names.arg=c("CS\nsingles", "pairs", "Bias\nsingles", "pairs"),
#         las=2,
#         ylim=c(0,600))


# window size is 7.99 x 5.22
# Save as PDF 8x4 inches


if(F){
  
  means = sapply(R, function(r)r[980,2])
  ers = (sapply(R, function(r)r[980,3])/2)^2
  best = which.min(means)
  pnorm( (means[best] - means)/sqrt(ers[best]+ers) ) >0.05
  
  means = sapply(R, function(r)r[480,2])
  ers = (sapply(R, function(r)r[480,3])/2)^2
  best = which.min(means)
  pnorm( (means[best] - means)/sqrt(ers[best]+ers) ) >0.05
  
}



