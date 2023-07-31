FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y wget curl unzip build-essential datamash graphviz

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.1.0-1-Linux-x86_64.sh \
	&& bash Miniconda3-py38_23.1.0-1-Linux-x86_64.sh -b -p /opt/conda \
	&& rm Miniconda3-py38_23.1.0-1-Linux-x86_64.sh
ENV PATH="${PATH}:/opt/conda/bin"

RUN mkdir -p /home/download
WORKDIR /home/download

RUN conda install -y -c bioconda seqkit minimap2
RUN apt-get install -y cmake libz-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev \
	libncursesw5-dev libcurl4-openssl-dev libffi-dev python3-all-dev libssl-dev python3-virtualenv

RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
RUN wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
RUN tar -xvf samtools-1.14.tar.bz2
RUN tar -xvf htslib-1.14.tar.bz2
RUN tar -xvf minimap2-2.17_x64-linux.tar.bz2
RUN tar -xvf bcftools-1.14.tar.bz2
RUN cd samtools-1.14 && make && cd ..
RUN cd htslib-1.14 && make && cd ..
RUN cd bcftools-1.14 && make && cd ..
ENV PATH="/home/download/samtools-1.14:/home/download/bcftools-1.14:/home/download/htslib-1.14:${PATH}"

RUN wget https://github.com/fenderglass/Flye/archive/refs/tags/2.9.1.tar.gz
RUN tar -xvf /home/download/2.9.1.tar.gz
RUN cd Flye-2.9.1/ && make && cd -
RUN ln -s /home/download/Flye-2.9.1/bin/flye /usr/local/bin/flye

RUN apt-get install -y git git-lfs python3.8-venv
RUN git clone https://github.com/nanoporetech/medaka.git
RUN cd /home/download/medaka && make install && cd -
ENV PATH="${PATH}:/home/download/medaka/venv/bin"

RUN wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz \
    && tar -xzf SPAdes-3.15.5-Linux.tar.gz

RUN apt-get install -y default-jre git
RUN wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
RUN git clone https://github.com/lh3/bwa.git \
	&& cd bwa && make \
	&& cp bwa /usr/local/bin \
	&& cd /home/download

RUN git clone https://github.com/rrwick/Unicycler.git \
	&& cd Unicycler && make \
	&& cd /home/download

RUN git clone https://github.com/lbcb-sci/racon \
	&& cd racon && mkdir build && cd build \
	&& cmake -DCMAKE_BUILD_TYPE=Release .. && make \
	&& cp /home/download/racon/build/bin/racon /usr/local/bin \
	&& cd /home/download

RUN git clone https://github.com/rrwick/Filtlong.git \
	&& cd Filtlong && make -j \
	&& cp bin/filtlong /usr/local/bin \
	&& cd /home/download

RUN apt-get install -y ncbi-blast+

RUN mkdir -p /home/scripts
WORKDIR /home/scripts
