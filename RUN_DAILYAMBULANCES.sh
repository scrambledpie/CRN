#!/bin/bash

Rscript DailyAmbulances_runner.R $(($1+1)) $2

rsync -rv EachData1/ huanan:/home/michael/OPUS/damb/$3
rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/damb/$3




