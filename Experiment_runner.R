# running within Rstudio
cat("Running locally \n")


# Optimization Run
method   = 5 # see list of optimizers below
problem  = 1 # the test function to use, see list below
BOseed   = 1 # the starting seed for the test problem
Ns0      = 5 # number of seeds to start the optimizer
Budget   = 25 # the total number of samples to take
filename = NULL # where to save results


# Make the TestFunction
source('TestFuns/TestFuns.R')
TESTFUNS = c(Build_Xie_ATO_cpp_Testfun,
             Build_Ambulance_Testfun,
             Build_DailyAmbulance_Testfun)

TestFun = TESTFUNS[[problem]](BOseed, numtestseeds=2000, runlength=1)[[1]]
ran = attr(TestFun, 'ran') # bounds of the search space.


# pick the optimizer from the list
source('CRN_BO/Optimizers.R')
source('CRN_BO/SOSA.R')
ALGORITHMS = c(BO_KG,                   # normal Knowledge Gradient (KG)
               BO_CRNKG_CS,             # KG-CRN with compound spheric noise model
               BO_PWKG_CS,              # KG-PW with compound spheric noise model
               BO_CRNKG_CSW,            # KG-CRN with full noise model: compound spheric + wiggles
               BO_PWKG_CSW,             # PW-KG with full noise model
               BO_CRNKG_CS_allseeds,    # KG-CRN with CS model, acq fun is optimized for each seed (slow)
               BO_CRNKG_CSW_allseeds,   # KG-CRN with full model, acq fun is optimized for each seed (slow)
               BO_CRNKG_CSW_allseeds)   # SOSA, the algrotihms proposed by 

Optim_class = ALGORITHMS[[method]]


# initialise the optimizer object
Optimizer = Optim_class$new(TestFun, ran, BOseed, myID=filename, rounding=T)


# exectute the optimizer
Optimizer$optimize(Budget0 = 20, Budget = Budget,
                   N0  = 1000, Na  = 10, maxevals  = 100,
                   PN0 = 4000, PNa = 20, Pmaxevals = 200)

cat("Finished and Saved ", filename)


  