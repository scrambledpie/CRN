setwd("/Users/pearce/CRN/")
source("plotting/plot_utils.R")

Results = readRDS("/Users/pearce/CopiedResults/8001")

M = sapply(Results, function(r)r$method)

R = Get_multiple_OC(Results, M)

S = Get_multiple_seed_reuse(Results, M)

HP_14 = Get_multiple_Hypers(Results, M, 14)
HP_15 = Get_multiple_Hypers(Results, M, 15)
HP_16 = Get_multiple_Hypers(Results, M, 16)

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

cols[1] = "blue" # KG
cols[2] = "red" # PWKG CS
cols[3] = "green" # faulty
cols[4] = "darkred" # PWKG+W
cols[5] = "darkgreen" # CRNKG PW
cols[6] = "green" #CRNKG CSW


par(mfrow=c(3,3))

# FIRST ROW
Res_plot(R, cols = cols, main = "Time to patient", erscale = 0)
Res_plot(R, cols = cols, main = "Time to patient")
Res_plot(S, cols = cols, main = "seed reuse")


# SECOND ROW
CRNbars = Get_seed_bars(Results[M==4], 30)
barplot(CRNbars, names.arg = 1:ncol(CRNbars), main="CRN CS", ylim=c(0, 200))

CRNbars = Get_seed_bars(Results[M==6], 30)
barplot(CRNbars, names.arg = 1:ncol(CRNbars), main="CRN wig", ylim=c(0, 200))

CS_PWbars = Get_seed_pairs_singles(Results[M==5])
CSW_PWbars = Get_seed_pairs_singles(Results[M==7])
barplot(c(CS_PWbars, CSW_PWbars), main="PW: CS, CSW", names.arg=c("CS singles", "CS pairs", "CSW singles", "CSW pairs"))
# barplot(CSW_PWbars, main="CSW_PW", names.arg=c("singles", "pairs"))

# THIRD ROW
Res_plot(HP_14, cols=cols, main="wiggle")
Res_plot(HP_15, cols=cols, main="offset")
Res_plot(HP_16, cols=cols, main="noise")





