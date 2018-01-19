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
        libdbd-mysql \
        libgsl0-dev \
        liblzma-dev \
        libmariadb-client-lgpl-dev \
        libncurses5-dev \
        libpcre3-dev \
        libreadline-dev \
        libssl-dev \
        libxml2-dev \
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
ENV TRIMGALORE_VERSION="0.4.5"
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/${TRIMGALORE_VERSION}.tar.gz -o /opt/trimgalore_${TRIMGALORE_VERSION}.tar.gz && \
    tar xvzf /opt/trimgalore_${TRIMGALORE_VERSION}.tar.gz -C /opt/ && \
    ln -s /opt/TrimGalore-${TRIMGALORE_VERSION}/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/trimgalore_${TRIMGALORE_VERSION}.tar.gz

# Install BWA
ENV BWA_VERSION="0.7.15"
RUN curl -fsSL https://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2 -o /opt/bwa-${BWA_VERSION}.tar.bz2 && \
    tar xvjf /opt/bwa-${BWA_VERSION}.tar.bz2 -C /opt/ && \
    cd /opt/bwa-${BWA_VERSION};make && \
    ln -s /opt/bwa-${BWA_VERSION}/bwa /usr/local/bin/bwa && \
    rm /opt/bwa-${BWA_VERSION}.tar.bz2

# Install SAMTools
ENV SAMTOOLS_VERSON="1.6"
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
ENV R_VERSION="R-3.4.2"
RUN curl -fsSL https://cran.r-project.org/src/base/R-3/${R_VERSION}.tar.gz -o /opt/${R_VERSION}.tar.gz && \
    tar xvzf /opt/${R_VERSION}.tar.gz -C /opt/ && \
    cd /opt/${R_VERSION};./configure;make;make install && \
    rm /opt/${R_VERSION}.tar.gz

# Install core R dependencies
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'https://ftp.acc.umu.se/mirror/CRAN/'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('caTools',dependencies=TRUE)" && \
    Rscript -e "install.packages('snow',dependencies=TRUE)" && \
    Rscript -e "install.packages('doMC',dependencies=TRUE)" && \
    Rscript -e "install.packages('utils',dependencies=TRUE)" && \
    Rscript -e "install.packages('stringr',dependencies=TRUE)" && \
    Rscript -e "install.packages('markdown',dependencies=TRUE)" && \
    Rscript -e "install.packages('evaluate',dependencies=TRUE)" && \
    Rscript -e "install.packages('ggplot2',dependencies=TRUE)" && \
    Rscript -e "install.packages('knitr',dependencies=TRUE)" && \
    Rscript -e "install.packages('RMySQL',dependencies=TRUE)"

# Install R Bioconductor packages
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("BSgenome", "Rsamtools", "ShortRead", "GenomicRanges", "GenomicFeatures", "ensembldb", "ChIPpeakAnno", "biomaRt", "rtracklayer", "BSgenome.Hsapiens.UCSC.hg19", "org.Hs.eg.db", "BSgenome.Mmusculus.UCSC.mm10", "org.Mm.eg.db"))' >> /opt/packages.r && \
    Rscript /opt/packages.r && \
    mkdir /usr/local/lib/R/site-library

# Install phantompeakqualtools
ENV SPP_VERSION="1.15"
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
ENV DEEPTOOLS_VERSION="2.5.4"
RUN pip install git+git://github.com/fidelram/deepTools.git@$DEEPTOOLS_VERSION

# Install ngsplot
ENV NGSPLOT_VERSION="2.63"
RUN curl -fsSL https://github.com/shenlab-sinai/ngsplot/archive/${NGSPLOT_VERSION}.tar.gz -o /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz && \
    tar xvzf /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz -C /opt/ && \
    rm /opt/ngsplot_${NGSPLOT_VERSION}.tar.gz
ENV PATH=${PATH}:/opt/ngsplot-${NGSPLOT_VERSION}/bin
ENV NGSPLOT=/opt/ngsplot-${NGSPLOT_VERSION}/
RUN curl -fsSL https://export.uppmax.uu.se/b2013064/test-data/ngi-chipseq_test_set.tar.bz2 -o /opt/ngi-chipseq_test_set.tar.bz2 && \
    tar xvjf /opt/ngi-chipseq_test_set.tar.bz2 -C /opt/ && \
    echo y | ngsplotdb.py install /opt/ngi-chipseq_test_set/ngsplotdb_mm10_75_3.00.tar.gz && \
    echo y | ngsplotdb.py install /opt/ngi-chipseq_test_set/ngsplotdb_hg19_75_3.00.tar.gz && \
    rm /opt/ngi-chipseq_test_set.tar.bz2 && \
    rm -rf /opt/ngi-chipseq_test_set

# Install MACS
RUN pip install MACS2

# Install MultiQC
ENV MULTIQC_VERSION="v1.4"
RUN pip install git+git://github.com/ewels/MultiQC.git@$MULTIQC_VERSION

# Create root directories for UPPMAX and c3se hebbe
RUN mkdir /pica /lupus /crex1 /crex2 /proj /scratch /sw \
          /c3se /local /apps
