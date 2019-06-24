
# Make the TestFunction
source('TestFuns/TestFuns.R')
TestFun = Build_Ambulance_Testfun(BOseed, 10000, 1, 10)[[1]]
ran = attr(TestFun, 'ran')


# pick the algortihm from the list
source('CRN_BO/Optimizers.R')


# exectute the optimizer
BB = lapply(1:40, function(BOseed){
  TestFun = Build_Ambulance_Testfun(BOseed, 10000, 1, 10)[[1]]
  ran = attr(TestFun, 'ran')
  AA = RANDOM_SEARCH$new(TestFun, ran)
  AA$optimize(Budget0 = 20, Budget = Budget, BOseed)
  
  return(AA)
})

saveRDS(BB, "Ambulances_random.rdata")
cat("Finished and Saved ", filename)


