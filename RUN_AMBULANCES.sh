#!/bin/bash

Rscript Ambulances_runner.R $1

rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/amb/$2

Rscript Ambulances_runner.R $(($1+1000))

rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/amb/$2


