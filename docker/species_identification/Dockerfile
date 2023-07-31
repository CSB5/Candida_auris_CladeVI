FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y wget curl unzip build-essential git libz-dev

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh \
	&& bash Miniconda3-py310_23.1.0-1-Linux-x86_64.sh -b -p /opt/conda \
	&& rm Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
ENV PATH="${PATH}:/opt/conda/bin"

RUN git clone https://github.com/Mykrobe-tools/mykrobe.git && cd mykrobe
RUN pip3 install . && mykrobe panels update_metadata && mykrobe panels update_species all

RUN mkdir -p /home/download
WORKDIR /home/download

RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
RUN tar -xvf v2.1.2.tar.gz
RUN cd kraken2-2.1.2 && ./install_kraken2.sh . && cd -

RUN mkdir -p /home/scripts
WORKDIR /home/scripts
