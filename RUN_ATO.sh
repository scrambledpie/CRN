#/bin/bash


Rscript ATO_runner.R $(($1 + 1))

rsync -rv EachData1/ jamon:/storage/maths/phrnaj/opus/ato/$2
rsync -rv EachData1/ huanan:/home/michael/OPUS/ato/$2


