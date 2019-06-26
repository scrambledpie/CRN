#!/bin/bash

Rscript DailyAmbulances_runner.R $(($1+1)) $2

rsync -rv EachData1/ huanan:/home/michael/OPUS/amb/$3
rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/amb/$3




