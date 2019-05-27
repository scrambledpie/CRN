
Args = commandArgs(trailingOnly = T)

BOseed = as.integer(Args[1]) + 1
Foldername = Args[2]

for(rho_i in 11:11){
  for(m in c(5,4,2)){
    cmd = paste("Rscript Vary_wiggle_runner.R", m, rho_i, BOseed)
    
    system(cmd)
    
  }
}

system(paste("rsync -rv EachData1/ huanan:/home/michael/OPUS/wig/",Foldername,"/", sep=""))
system(paste("rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/wig/",Foldername,"/",sep=""))