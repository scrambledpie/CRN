if (CPU=="huanan")        setwd("~/Dropbox/PhD/CRN/")
if (grepl("Book", CPU))   setwd("/Volumes/DataStorage/Dropbox/PhD/CRN/git_CRN/")
if (grepl("pearce", CPU)) setwd("/Users/pearce/CRN/")

source('TestFuns/TestFuns.R')

TestFuns = Build_DailyAmbulance_Testfun(BOseed, numtestseeds=10000, runlength=1, Simtime=1800)


par(mfrow=c(4,4))
for (i in 1:16){
  x1 = 0.1 + 0.8*(1 + c(cos(pi*i/16), sin(pi*i/16)) )/2
  x2 = 0.1 + 0.8*(1 + c(-cos(pi*i/16), -sin(pi*i/16)))/2
  print(x1, x2)
  ambulance1 = function(x, s=2){
    TestFuns[[1]](matrix(c(x, x1, x2, s), 1))  
  }
  ambulance1(c(0.5,0.5))
  
  x = seq(0, 1, len=20)
  X = cbind(x, rep(x, each=length(x)))
  
  
  Y = apply(X, 1, ambulance1)
  
  # par(mfrow=c(1,1))
  contour(x, x, matrix(Y, 20, 20))
  points(c(x1[1], x2[1]), c(x1[2], x2[2]), pch=19, cex=2)
  
  Streams = TestFuns[[4]](2)
  Calls = Streams$CallTimes[[1]]<1800
  CallLocs = Streams$CallLocs[[1]][Calls,]
  points(CallLocs)
  
  # Sys.sleep(0.5)
}