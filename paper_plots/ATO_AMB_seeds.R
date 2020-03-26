
# First load up the test functions
source("/Users/academic/CRN/TestFuns/TestFuns.R")
setwd("/Users/academic/CRN/")
par(mfrow=c(1,1))

# Make the TestFunction
source('TestFuns/TestFuns.R')
TESTFUNS = c(Build_Xie_ATO_cpp_Testfun,
             Build_Ambulance_Testfun,
             Build_DailyAmbulance_Testfun)

BOseed=1

problem = 1
TestFun = TESTFUNS[[problem]](BOseed, numtestseeds=2000, runlength=1)[[1]]
ran = attr(TestFun, 'ran') # bounds of the search space.

# make the ATO function
ATO_uniform = function(x, s) TestFun(c(rep(x, 8), s))
ATO_output = sapply(1:100, function(s){
  sapply(0:20, function(x){
    ATO_uniform(x,s)
  })
})

# initialize plot and plot each seed
plot(c(0, 20), range(ATO_output), col="white", main="ATO")
for(s in 1:6){
 lines(0:20, ATO_output[, s], col="grey") 
 points(0:20, ATO_output[, s], pch=19, col="grey") 
}
lines(0:20, apply(ATO_output, 1, mean), lwd=3)


# make the Ambulances function
problem = 3
AMBFun = TESTFUNS[[problem]](BOseed, numtestseeds=2000, runlength=1)[[1]]
AMB_diag = function(x, s) AMBFun(c(18, 8, x[1], 
                                   10, 15, x[1], s))
AMB_output = sapply(1:10, function(s){
  sapply((0:100)/5, function(x){
    abs(AMB_diag(x,s))
  })
})

# initialize plotand plot each seed
plot(c(0, 20), c(0.125, 0.19), col="white", main="AIS")
for(s in 1:4){
  lines((0:100)/5, AMB_output[, s], col="grey")
}
lines((0:100)/5, apply(AMB_output, 1, mean), lwd=3)

