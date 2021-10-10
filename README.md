# Bayesian Optimization Allowing for Common Random Numbers

The code for producucing the result in the papaer to appear in the journal "Operations Research". For further details check out
- Blog post [BO-CRN Blog Post](https://bayesianblog.com/BO-CRN/) 
- arxiv preprint [Bayesian Optimization Allowing for Common Random Numbers](https://arxiv.org/abs/1910.09259)

A few disclaimers
- the code is written entirely in R (yes, python would be a loot better, hopefully one day!)
- the code was written by a PhD student (me in 2018) with very little knowledge of good coding practice at the time!

## Running code

Clone the repo and step inside
```
git clone https://github.com/scrambledpie/CRN.git
cd CRN
```
To run different baselines and differnt methods, simply edit the settings in the `Experiment_runner.R` file.

### Run from R terminal
Assuming you have a version or R installed on your computer.
1. Install R packages `FastGP`, `Rcpp`, `R6`, `testit`
```
install.packages(c("FastGP", "Rcpp", "R6", "testit"), repos="https://cloud.r-project.org")
```
2. Run an example experiment from the command line
```
source('Experiment_runner.R')
```

### Run from bash
Assuming you have a version or R installed on your computer.
1. Install R packages `FastGP`, `Rcpp`, `R6`, `testit`
```
Rscript -e 'install.packages(c("FastGP", "Rcpp", "R6", "testit"), repos="https://cloud.r-project.org")'
```
2. Run an example experiment from the command line
```
Rscript Experiment_runner.R
```

### Docker
If you do not have/want R installed, then use the docker container.
1. Build the docker image form the Dockerfile
```
docker build -f Dockerfile -t bo_crn .
```
2. Mount the local code dir into a container and run it
```
docker run -v $(pwd):/code/ bo_crn Rscript Experiment_runner.R
```


