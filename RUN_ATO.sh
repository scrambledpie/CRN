#/bin/bash

Rscript ATO_runner.R $(($1 + 1001))

rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/ato/$2

Rscript ATO_runner.R $(($1 + 1))

rsync -rv EachData1/ godzilla:/storage/maths/phrnaj/opus/ato/$2


