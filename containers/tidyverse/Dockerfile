FROM mambaorg/micromamba:1.5.8

LABEL software="R with Bioconductor KEGGREST"
LABEL version="4.3.2"

USER root

RUN apt-get update && apt-get install -y \
    procps \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libblas-dev liblapack-dev gfortran \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Installer R et packages via micromamba dans l'environnement base
RUN micromamba install -y -n base -c conda-forge -c bioconda -c r \
    r-base=4.3.2 \
    bioconductor-keggrest \
    r-tidyverse \
    r-venndiagram \
    r-ggrepel \
    r-scales \
    && micromamba clean -a -y

ENV PATH="/opt/conda/envs/base/bin:/opt/conda/bin:$PATH"

WORKDIR /usr/src/app

CMD ["R"]
