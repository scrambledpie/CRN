FROM r-base:latest

RUN Rscript -e 'install.packages("FastGP", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("Rcpp", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages(c("R6", "testit"), repos="https://cloud.r-project.org")'

WORKDIR /code

# ENTRYPOINT ["Rscript"]
# CMD ["-e", "'print(\"give a script\")'"]