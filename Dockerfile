FROM nfcore/base
MAINTAINER Phil Ewels <phil.ewels@scilifelab.se>
  LABEL authors="phil.ewels@scilifelab.se" \
description="Docker image containing all requirements for the nfcore/ChIPseq pipeline"

COPY environment.yml /
  RUN conda update -n base conda
RUN conda update --all
RUN conda env update -n root --file environment.yml  && conda clean -a
ENV PATH /opt/conda/envs/nfcore-chipseq-1.4dev/bin:$PATH

ENV NGSPLOT_VERSION="2.63"
RUN curl -fsSL https://github.com/shenlab-sinai/ngsplot/archive/${NGSPLOT_VERSION}.tar.gz -o /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz && \
    tar xvzf /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz -C /opt/ && \
    rm /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz
ENV PATH=${PATH}:/opt/ngsplot-${NGSPLOT_VERSION}/bin
ENV NGSPLOT=/opt/ngsplot-${NGSPLOT_VERSION}/

RUN wget "https://drive.google.com/uc?export=download&id=0B5hDZ2BucCI6SURYWW5XdUxnbW8" -O ngsplotdb_hg19_75_3.00.tar.gz && \
    echo y | ngsplotdb.py install ngsplotdb_hg19_75_3.00.tar.gz && \
    rm -rf  ngsplotdb_hg19_75_3.00.tar.gz && \
    wget "https://drive.google.com/uc?export=download&id=0B5hDZ2BucCI6S3E4dVprdlF2YW8" -O ngsplotdb_hg38_76_3.00.tar.gz && \
    echo y | ngsplotdb.py install ngsplotdb_hg38_76_3.00.tar.gz && \
    rm -rf  ngsplotdb_hg38_76_3.00.tar.gz && \
    wget "https://drive.google.com/uc?export=download&id=0B5hDZ2BucCI6NXNzNjZveXdadU0" -O ngsplotdb_mm10_75_3.00.tar.gz && \
    echo y | ngsplotdb.py install ngsplotdb_mm10_75_3.00.tar.gz && \
    rm -rf  ngsplotdb_mm10_75_3.00.tar.gz 

RUN git clone https://github.com/kundajelab/phantompeakqualtools  && \
    mv phantompeakqualtools /opt/  && \
    echo 'library(caTools)' | cat - /opt/phantompeakqualtools/run_spp.R > temp && mv temp /opt/phantompeakqualtools/run_spp.R && \
    chmod 755 /opt/phantompeakqualtools/* && \
    echo 'alias run_spp.R="Rscript /opt/phantompeakqualtools/run_spp.R"' >> ~/.bashrc 
ENV PATH=${PATH}:/opt/phantompeakqualtools

