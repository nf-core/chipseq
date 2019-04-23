FROM nfcore/base
LABEL authors="phil.ewels@scilifelab.se" \
      description="Docker image containing all requirements for the nfcore/chipseq pipeline"

## GAWK has the 'and' function, needed for run_spp.R script from phantompeakqualtools package
RUN apt-get update && apt-get install -y gawk && apt-get clean -y
RUN echo 'alias awk=gawk' >> ~/.bashrc

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-chipseq-1.0dev/bin:$PATH
