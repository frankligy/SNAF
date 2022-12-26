#!/bin/bash


STAR --runMode genomeGenerate \
     --genomeDir /tmp/GenomeRef \
     --genomeFastaFiles /tmp/GRCh38.p13.genome.fa \
     --runThreadN 8 \
     --sjdbGTFfile /tmp/gencode.v36.annotation.gtf \
     --sjdbOverhang 100