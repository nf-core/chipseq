FROM nfcore/base:1.13.1
LABEL authors="Philip Ewels" \
      description="Docker image containing all software requirements for the nf-core/chipseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-chipseq-1.3.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-chipseq-1.3.0dev > nf-core-chipseq-1.3.0dev.yml
