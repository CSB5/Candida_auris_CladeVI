FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]

RUN apt-get update
RUN apt-get install -y wget curl unzip build-essential 

# install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.1.0-1-Linux-x86_64.sh \
	&& bash Miniconda3-py38_23.1.0-1-Linux-x86_64.sh -b -p /opt/conda \
	&& rm Miniconda3-py38_23.1.0-1-Linux-x86_64.sh
ENV PATH="${PATH}:/opt/conda/bin"

# Install mash
RUN mkdir -p /home/download
WORKDIR /home/download
RUN wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar \
    && tar -xvf mash-Linux64-v2.3.tar \
    && cp mash-Linux64-v2.3/mash /usr/local/bin/ \
    && rm -r mash-Linux64-v2.3 \
    && rm mash-Linux64-v2.3.tar

# Install snpEff and build C auris annotation database
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y default-jre-headless
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
    && unzip snpEff_latest_core.zip
RUN mkdir /home/download/snpEff/data \
    && cd /home/download/snpEff/data \
    && curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_002759435.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_002759435.2.zip" -H "Accept: application/zip" \
    && unzip GCA_002759435.2.zip \
    && mkdir candida_auris_GCA_002759435.2 \
    && mv ncbi_dataset/data/GCA_002759435.2/* candida_auris_GCA_002759435.2/ \
    && rm -r ncbi_dataset/ \
    && cd candida_auris_GCA_002759435.2/ \
    && ln -s genomic.gff genes.gff \
    && ln -s GCA_002759435.2_Cand_auris_B8441_V2_genomic.fna sequences.fa \
    && echo -e "\ncandida_auris_GCA_002759435.2.genome : candida_auris_GCA_002759435.2\n" >> /home/download/snpEff/snpEff.config \
    && java -jar /home/download/snpEff/snpEff.jar build -gff3 -v candida_auris_GCA_002759435.2 -noCheckCds -noCheckProtein

WORKDIR /home
