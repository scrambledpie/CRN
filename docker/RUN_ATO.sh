#/bin/bash


Rscript ATO_runner.R $(($1 + 1)) $2

rsync -rv EachData1/ huanan:/home/michael/OPUS/ato/$3
rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/ato/$3



