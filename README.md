# WGS analysis for *Candida auris*
This repository contains the source code of the Nextflow pipelines for sequencing analyses, as well as a Jupyter Notebook implementing a machine-learning model for *Candida auris* clade detection.

## Prerequisite
The dependencies required before running the workflow are following
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#requirements)
- [Docker](https://docs.docker.com/engine/install/)
- All the software used in this paper was packaged inside the docker image. Please refer to `docker` directory to build the necessary containers for the workflow

## Workflow description
* `candida_auris_isolate_illumina.nf`<br>
    The workflow processes **Illumina reads**, following these steps<br>
    1. Preprocess the reads by trimming the adapters and quality filtering using [Fastp](https://github.com/OpenGene/fastp)
    2. Calculate statistics information using [Fastqc](https://github.com/s-andrews/FastQC) and [Seqkit](https://github.com/shenwei356/seqkit)
    3. Perform reads-based species identification using [Mash](https://github.com/marbl/Mash)
    4. Perform sequencing depth checking using our custom script
    5. Perform genome assembly and filtering short contigs using [SPAdes](https://github.com/ablab/spades) and [Seqkit](https://github.com/shenwei356/seqkit)
    6. Perform assembly QC evaluation using [Quast](https://github.com/ablab/quast)
    7. Perform contig species identification using [Mash](https://github.com/marbl/Mash)
    8. Perform variant calling using [SnpEff](https://github.com/pcingola/SnpEff) and [Snippy](https://github.com/tseemann/snippy)
- `candida_auris_isolate_nanopore.nf`<br>
    The workflow processes a **Nanopore read**, following these steps<br>
    1. Preprocess the reads by quality filtering using [chopper](https://github.com/wdecoster/chopper)
    2. Calculate statistics information using [NanoStat](https://github.com/wdecoster/nanostat) and [Seqkit](https://github.com/shenwei356/seqkit)
    3. Perform sequencing depth checking using our custom script
    4. Perform genome assembly and draft consensus genome construction using [Flye](https://github.com/fenderglass/Flye)
    5. Perform genome polishing and short contigs filtering using [Medaka](https://github.com/nanoporetech/medaka) and [Seqkit](https://github.com/shenwei356/seqkit)
    6. Perform assembly QC evaluation using [Quast](https://github.com/ablab/quast) with `--nanopore` option
- `candida_auris_isolate_illumina_polish.nf`<br>
    The workflow performs genome polishing by using Illumina short reads to improve the quality of the nanopore-based assembled genome, following these steps.<br>
    1. Perform long-read indexing, then align short-reads to it using [BWA](https://github.com/lh3/bwa)
    2. From i. we generate a sorted index file in BAM format using [SAMtools](https://github.com/samtools/samtools)
    3. Perform genome polishing from the long-read and its index file using [Pilon](https://github.com/broadinstitute/pilon)
    4. Perform assembly QC evaluation using [Quast](https://github.com/ablab/quast) with `--nanopore` option

## Running the workflow
The below code is an example command to run the workflow<br><br>
`candida_auris_isolate_illumina.nf`<br>
Note that the Illumina paired-end read was retrieved from the SRA database by using `--sra_accession` option customized in our script.<br> For the workflow configuration, please refer to `nextflow.config`
```bash
# =================================
# candida_auris_isolate_illumina.nf
# The Nextflow script for Illumina
# =================================

nextflow candida_auris_isolate_illumina.nf \
-with-report /data/sgh_candida_auris/20230630_illumina_analysis/SRR24877249/SRR24877249.nextflow.report.html \
-with-trace /data/sgh_candida_auris/20230630_illumina_analysis/SRR24877249/SRR24877249.nextflow.trace.txt \
-with-timeline /data/sgh_candida_auris/20230630_illumina_analysis/SRR24877249/SRR24877249.nextflow.timeline.html \
-w /data/sgh_candida_auris/20230630_illumina_analysis/nextflow_work \
--outdir /data/sgh_candida_auris/20230630_illumina_analysis/SRR24877249 \
--sample_id SRR24877249 \
--genome_size 12500000 \
--sra_accession SRR24877249
```

`candida_auris_isolate_nanopore.nf`<br>
For the workflow configuration, please refer to `nextflow.config`
```bash
# =================================
# candida_auris_isolate_nanopore.nf
# The Nextflow script for Nanopore
# =================================

nextflow candida_auris_isolate_nanopore.nf \
-with-report /home/ubuntu/data/sgh_candida_auris/20230604_analysis_results/N00521_BC17/nextflow_report.html \
-with-trace /home/ubuntu/data/sgh_candida_auris/20230604_analysis_results/N00521_BC17/nextflow_trace.txt \
-with-timeline /home/ubuntu/data/sgh_candida_auris/20230604_analysis_results/N00521_BC17/nextflow_timeline.html \
-w /home/ubuntu/data/sgh_candida_auris/20230604_analysis_results/nextflow_work \
--outdir /home/ubuntu/data/sgh_candida_auris/20230604_analysis_results/N00521_BC17 \
--sample_id N00521_BC17 \
--nanopore_reads /home/ubuntu/data/sgh_candida_auris/N00521_BC17.fastq.gz 
```

`candida_auris_isolate_illumina_polish.nf`<br>
For the workflow configuration, please refer to `nextflow.config`
```bash
# =================================
# candida_auris_isolate_illumina_polish.nf
# The Nextflow script for Hybrid genome assembly
# =================================

nextflow candida_auris_isolate_illumina_polish.nf \
-with-report /data/sgh_candida_auris/20230619_illumina_polish/F01567/nextflow_report.html \
-with-trace /data/sgh_candida_auris/20230619_illumina_polish/F01567/nextflow_trace.txt \
-with-timeline /data/sgh_candida_auris/20230619_illumina_polish/F01567/nextflow_timeline.html \
-w /data/sgh_candida_auris/20230619_illumina_polish/nextflow_work \
--outdir /data/sgh_candida_auris/20230619_illumina_polish/F01567 \
--sample_id F01567 \
--genome_size 12500000 \
--nanopore_reads /data/sgh_candida_auris/20230604_analysis_results/N00466_BC05/N00466_BC05_for_downstream.fastq \
--nanopore_contigs /data/sgh_candida_auris/20230604_analysis_results/N00466_BC05/N00466_BC05.nanopore.flye.medaka_x2.fasta \
--illumina_reads_1 /data/sgh_candida_auris/illumina_fastq/WMB1897_DKDL220004524-1a-AK17215-AK4966_HHJJGCCX2_L1_1.fq.gz \
--illumina_reads_2 /data/sgh_candida_auris/illumina_fastq/WMB1897_DKDL220004524-1a-AK17215-AK4966_HHJJGCCX2_L1_2.fq.gz 
```
## Data availability
Illumina and Nanopore sequencing data of three Clade VI isolates have been deposited in the National Centre for Biotechnology Information (NCBI) [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/docs/) under BioProject accession number [PRJNA1000034](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1000034).
<br><br>
Illumina data for large cohort analysis, a total of 4,475 unique NCBI accession ids were retrieved from the SRA database [The Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/docs/).<br>

## Machine-learning models
As a *proof-of-concept*, we provide a Jupyter notebook `machine_learning_for_clade_detection.ipynb` for the automatic detection of a new *Candida auris*. These approaches can potentially enhance genomic surveillance by the early identification and investigation of outlier genomes.
<br><br>
In short, Bayesian logistic regression models were trained based on SNP distances and previously reported clade information to learn a threshold for predicting whether a pair of genomes are from the same clade. The threshold was then used to determine the relationships between genome pairs (edge) in a graph, which captured clusters (connected components) representing existing and potential new clades. The analysis was based on 3,651 publicly available WGS and three Clade VI ([PRJNA1000034](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1000034)), of which 1,132 (31%) had previously reported clade information were included (A large SNP distance matrix file is available upon request.) At each time point, a graph was generated, where nodes represent genomes and edges link between two genomes that were predicted to belong to the same clade. The number of clusters (connected components) present in the graph represents the total number of clades predicted to be present in the dataset. An overview of our machine learning approach for detecting potential _Candida auris_ new clade is depicted in the figure below.

<img width="798" alt="image" src="https://github.com/kitanu-sg/discovery_of_the_6th_Candida_auris_clade/assets/11409990/2028297f-bee8-4efb-9719-096de83f8149">

