#!/bin/bash

Rscript Ambulances_runner.R $(($1+1))

rsync -rv EachData1/ jamon:/storage/maths/phrnaj/opus/amb/$2
rsync -rv EachData1/ huanan:/home/michael/OPUS/amb/$2


Rscript Ambulances_runner.R $(($1+801))

rsync -rv EachData1/ jamon:/storage/maths/phrnaj/opus/amb/$2
rsync -rv EachData1/ huanan:/home/michael/OPUS/amb/$2


