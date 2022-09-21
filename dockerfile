FROM quay.io/singlecellpipelinetest/miniconda3:4.8.2


RUN apt-get update --allow-releaseinfo-change && apt install build-essential gcc libfontconfig1 samtools -y && rm -rf /var/lib/apt/lists/*

RUN conda install -c bioconda minimap picard sambamba

RUN pip install git+https://github.com/diljotgrewal/ont_alignment:master
