FROM openjdk:8

LABEL author="Phil Ewels" \
    description="Docker image containing all requirements for NGI-ChIPseq pipeline" \
    maintainer="phil.ewels@scilifelab.se"

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        g++ \
        gawk \
        gcc \
        gfortran \
        libboost-all-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgsl0-dev \
        liblzma-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        make \
        python-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Use gawk instead of awk (needed for phantompeakqualtools)
RUN update-alternatives --set awk /usr/bin/gawk

# Install pip
RUN curl -fsSL https://bootstrap.pypa.io/get-pip.py -o /opt/get-pip.py && \
    python /opt/get-pip.py && \
    rm /opt/get-pip.py

# Install FastQC
ENV FASTQC_BIN="fastqc_v0.11.5.zip"
RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$FASTQC_BIN -o /opt/$FASTQC_BIN && \
    unzip /opt/$FASTQC_BIN -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/$FASTQC_BIN

# Install cutadapt
RUN pip install cutadapt

# Install TrimGalore
ENV TRIMGALORE_BIN="trim_galore_v0.4.2.zip"
RUN mkdir /opt/TrimGalore && \
    curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/$TRIMGALORE_BIN -o /opt/TrimGalore/$TRIMGALORE_BIN && \
    unzip /opt/TrimGalore/$TRIMGALORE_BIN -d /opt/TrimGalore && \
    ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/TrimGalore/$TRIMGALORE_BIN

# Install BWA
ENV BWA_VERSION="0.7.15"
RUN curl -fsSL https://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2 -o /opt/bwa-${BWA_VERSION}.tar.bz2 && \
    tar xvjf /opt/bwa-${BWA_VERSION}.tar.bz2 -C /opt/ && \
    cd /opt/bwa-${BWA_VERSION};make && \
    ln -s /opt/bwa-${BWA_VERSION}/bwa /usr/local/bin/bwa && \
    rm /opt/bwa-${BWA_VERSION}.tar.bz2

# Install SAMTools
ENV SAMTOOLS_VERSON="1.4"
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSON}/samtools-${SAMTOOLS_VERSON}.tar.bz2 -o /opt/samtools-${SAMTOOLS_VERSON}.tar.bz2 && \
    tar xvjf /opt/samtools-${SAMTOOLS_VERSON}.tar.bz2 -C /opt/ && \
    cd /opt/samtools-${SAMTOOLS_VERSON};make;make install && \
    rm /opt/samtools-${SAMTOOLS_VERSON}.tar.bz2

# Install PicardTools
ENV PICARD_VERSION="picard-tools-2.0.1"
RUN curl -fsSL https://github.com/broadinstitute/picard/releases/download/2.0.1/${PICARD_VERSION}.zip -o /opt/${PICARD_VERSION}.zip && \
    unzip /opt/${PICARD_VERSION}.zip -d /opt/ && \
    rm /opt/${PICARD_VERSION}.zip
ENV PICARD_HOME /opt/${PICARD_VERSION}

# Install BEDTools
ENV BEDTOOLS_VERSION="bedtools-2.26.0"
RUN curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v2.26.0/${BEDTOOLS_VERSION}.tar.gz -o /opt/${BEDTOOLS_VERSION}.tar.gz && \
    tar xvzf /opt/${BEDTOOLS_VERSION}.tar.gz -C /opt/ && \
    cd /opt/bedtools2; make && \
    cp /opt/bedtools2/bin/* /usr/local/bin/ && \
    rm /opt/${BEDTOOLS_VERSION}.tar.gz

# Install R
ENV R_VERSION="R-3.3.3"
RUN curl -fsSL https://cran.r-project.org/src/base/R-3/${R_VERSION}.tar.gz -o /opt/${R_VERSION}.tar.gz && \
    tar xvzf /opt/${R_VERSION}.tar.gz -C /opt/ && \
    cd /opt/${R_VERSION};./configure;make;make install && \
    rm /opt/${R_VERSION}.tar.gz

# Install core R dependencies
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://ftp.acc.umu.se/mirror/CRAN/'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('caTools',dependencies=TRUE)" && \
    Rscript -e "install.packages('snow',dependencies=TRUE)" && \
    Rscript -e "install.packages('doMC',dependencies=TRUE)" && \
    Rscript -e "install.packages('utils',dependencies=TRUE)"

# Install R Bioconductor packages
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("BSgenome", "Rsamtools", "ShortRead"))' >> /opt/packages.r && \
    Rscript /opt/packages.r && \
    mkdir /usr/local/lib/R/site-library

# Install phantompeakqualtools
ENV SPP_VERSION="1.14"
ENV PHANTOMPEAKQUALTOOLS_VERSION="v.1.1"
RUN curl -fsSL https://github.com/hms-dbmi/spp/archive/${SPP_VERSION}.tar.gz -o /opt/SPP_${SPP_VERSION}.tar.gz && \
    Rscript -e "install.packages('/opt/SPP_${SPP_VERSION}.tar.gz',dependencies=TRUE)" && \
    rm /opt/SPP_${SPP_VERSION}.tar.gz && \
    curl -fsSL https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.${PHANTOMPEAKQUALTOOLS_VERSION}.tar.gz -o /opt/phantompeakqualtools.${PHANTOMPEAKQUALTOOLS_VERSION}.tar.gz && \
    tar xvzf /opt/phantompeakqualtools.${PHANTOMPEAKQUALTOOLS_VERSION}.tar.gz -C /opt/ && \
    chmod 755 /opt/phantompeakqualtools/* && \
    echo 'alias run_spp.R="Rscript /opt/phantompeakqualtools/run_spp.R"' >> ~/.bashrc && \
    rm /opt/phantompeakqualtools.${PHANTOMPEAKQUALTOOLS_VERSION}.tar.gz
ENV PATH=${PATH}:/opt/phantompeakqualtools

# Install DeepTools
ENV DEEPTOOLS_VERSION="2.4.3"
RUN pip install git+git://github.com/fidelram/deepTools.git@$DEEPTOOLS_VERSION

# Install ngsplot
ENV NGSPLOT_VERSION="2.61"
RUN curl -fsSL https://github.com/shenlab-sinai/ngsplot/archive/${NGSPLOT_VERSION}.tar.gz -o /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz && \
    tar xvzf /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz -C /opt/ && \
    rm /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz
ENV PATH=${PATH}:/opt/ngsplot-${NGSPLOT_VERSION}/bin
ENV NGSPLOT=/opt/ngsplot-${NGSPLOT_VERSION}/

# Install MACS
RUN pip install MACS2

# Install MultiQC
# ENV MULTIQC_VERSION v0.9
ENV MULTIQC_VERSION master
RUN pip install git+git://github.com/ewels/MultiQC.git@$MULTIQC_VERSION
