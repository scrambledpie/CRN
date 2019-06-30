rm(list=ls())
CPU = system("hostname",intern=T)
Args = commandArgs(trailingOnly = T)
debug = length(Args)==0
.libPaths(c(.libPaths(),"~/R/"))


if (CPU=="huanan")        setwd("~/Dropbox/PhD/CRN/")
if (grepl("Book", CPU))   setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/git_CRN/")
if (grepl("pearce", CPU)) setwd("/Users/pearce/CRN/")
  
  # Optimization Run
  method   = 7
  BOseed   = 1


# Make the TestFunction
source('TestFuns/TestFuns.R')
# TestFun = Build_Ambulance_Testfun(BOseed, numtestseeds=60000, runlength=1, NumCalls=5)[[1]]
TestFun = Build_DailyAmbulance_Testfun(BOseed, numtestseeds=60000, runlength=1, Simtime=1800)[[1]]
ran = attr(TestFun, 'ran')



cat(TestFun(c(0.2,0.2, 0.5, 0.5, 0.8, 0.8, 0) ))

