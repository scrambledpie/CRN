
Args = commandArgs(trailingOnly = T)

jobID = as.integer(Args[1]) + 1

BOseeds = rep(1:800, 2)
BOseed = BOseeds[jobID]
Foldername = Args[2]

if(jobID<=800){
  
  for(rho_i in c(1:5)){
    for(m in c(5,4,2)){
      cmd = paste("Rscript Vary_wiggle_runner.R", m, rho_i, BOseed)
      system(cmd)
    }
  }

}else{
  
  for(rho_i in c(6:11)){
    for(m in c(5,4,2)){
      cmd = paste("Rscript Vary_wiggle_runner.R", m, rho_i, BOseed)
      system(cmd)
    }
  }
  
}

system(paste("rsync -rv EachData1/ huanan:/home/michael/OPUS/wig/",Foldername,"/", sep=""))
system(paste("rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/wig/",Foldername,"/",sep=""))
system(paste("rsync -rv EachData1/ olcay:/home/pearce/OPUS/wig/",Foldername,"/", sep=""))