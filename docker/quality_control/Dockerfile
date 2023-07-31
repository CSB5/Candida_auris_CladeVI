FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y wget curl unzip build-essential datamash graphviz \
    cmake libz-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev \
    libncursesw5-dev libcurl4-openssl-dev libffi-dev python3-all-dev libssl-dev

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh \
    && bash Miniconda3-py39_4.11.0-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-py39_4.11.0-Linux-x86_64.sh
ENV PATH="${PATH}:/opt/conda/bin"

RUN apt-get install -y default-jre
RUN conda install -y -c bioconda seqkit bcftools minimap2
RUN pip install pandas networkx taxonomy ncbi-datasets-pylib scipy biopython tabulate nanostat

# install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip fastqc_v0.12.1.zip \
    && ln -s /FastQC/fastqc /usr/local/bin/fastqc

# install fastp
RUN wget http://opengene.org/fastp/fastp.0.23.2 \
    && mv fastp.0.23.2 fastp \
    && chmod a+x ./fastp \
    && mv /fastp /usr/local/bin

# install Bowtie 2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.0/bowtie2-2.5.0-linux-x86_64.zip \
    && unzip bowtie2-2.5.0-linux-x86_64.zip \
    && cp bowtie2-2.5.0-linux-x86_64/bowtie2* /usr/local/bin \
    && rm bowtie2-2.5.0-linux-x86_64.zip

# install minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    && tar -jxvf minimap2-2.24_x64-linux.tar.bz2 \
    && cp minimap2-2.24_x64-linux/* /usr/local/bin \
    && rm -r minimap2-2.24_x64-linux \
    && rm minimap2-2.24_x64-linux.tar.bz2

RUN mkdir -p /home/download
WORKDIR /home/download
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
RUN wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
RUN tar -xvf samtools-1.14.tar.bz2
RUN tar -xvf htslib-1.14.tar.bz2
RUN tar -xvf bcftools-1.14.tar.bz2
RUN cd samtools-1.14 && make && cd ..
RUN cd htslib-1.14 && make && cd ..
RUN cd bcftools-1.14 && make && cd ..
ENV PATH="/home/download/samtools-1.14:/home/download/bcftools-1.14:/home/download/htslib-1.14:${PATH}"

RUN mkdir -p /home/scripts
WORKDIR /home/scripts

RUN pip install pysam

COPY *.py /home/scripts/
