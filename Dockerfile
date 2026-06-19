FROM bioconductor/bioconductor_docker:RELEASE_3_18

WORKDIR /rnaseq

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "BiocManager::install(c('DESeq2', 'clusterProfiler', 'org.Mm.eg.db', 'tximport'), ask = FALSE, update = FALSE)"

RUN R -e "install.packages('ggplot2', repos='https://cloud.r-project.org')"

COPY R/ ./R/
COPY analysis.R ./

CMD ["Rscript", "-e", "source('analysis.R')"]`
