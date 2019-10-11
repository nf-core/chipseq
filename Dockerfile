FROM nfcore/base:1.7
LABEL authors="Philip Ewels" \
      description="Docker image containing all requirements for nf-core/chipseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-chipseq-1.0.1dev/bin:$PATH
