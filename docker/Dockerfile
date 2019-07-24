FROM rocker/verse:3.6.1

RUN apt-get update
RUN apt-get -y upgrade

RUN apt-get -y install libudunits2-dev
RUN apt-get -y install graphviz

RUN Rscript -e "update.packages(ask=FALSE)"
RUN Rscript -e "install.packages('reshape2', quiet = TRUE)"
RUN Rscript -e "install.packages('stringapi', quiet = TRUE)"
RUN Rscript -e "install.packages('tidygraph', quiet = TRUE)"
RUN Rscript -e "install.packages('tidyverse', quiet = TRUE)"
RUN Rscript -e "install.packages('xtable', quiet = TRUE)"
RUN Rscript -e "install.packages('GGally', quiet = TRUE)"
RUN Rscript -e "install.packages('ggraph', quiet = TRUE)"
RUN Rscript -e "install.packages('gghighlight', quiet = TRUE)"
RUN Rscript -e "install.packages('cowplot', quiet = TRUE)"
RUN Rscript -e "install.packages('circlize', quiet = TRUE)"
RUN Rscript -e "install.packages('BiocManager', quiet = TRUE)"

RUN Rscript -e "BiocManager::install('ChIPseeker', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('ComplexHeatmap', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('DESeq2', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('DGCA', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('EnrichedHeatmap', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('ExperimentHub', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('GenomicRanges', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('Gviz', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('SummarizedExperiment', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('TCseq', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('curatedAdipoRNA', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('org.Mm.eg.db', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('rtracklayer', quiet = TRUE)"
RUN Rscript -e "BiocManager::install('curatedAdipoChIP', quiet = TRUE)"
