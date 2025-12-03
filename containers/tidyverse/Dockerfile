FROM rocker/r-base:4.3.2

RUN apt-get update && apt-get install -y \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libblas-dev liblapack-dev gfortran \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('tidyverse','VennDiagram','ggrepel','BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('KEGGREST','DESeq2'), update=FALSE, ask=FALSE)"

WORKDIR /usr/src/app

CMD ["R"]