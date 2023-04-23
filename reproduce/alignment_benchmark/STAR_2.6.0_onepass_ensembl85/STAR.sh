#!/bin/bash

FASTQ1=$1
FASTQ2=${FASTQ1/read1/read2}
SAMPLE=$(basename $FASTQ1 .fq.gz)

DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -M 128000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR
module load STAR/2.6.1

STAR --genomeDir /data/salomonis2/Genomes/Star-Index-GRCH38/Grch38-STAR-index --readFilesIn $FASTQ1 $FASTQ2 --readFilesCommand gunzip -c --outFileNamePrefix $DIR/$SAMPLE --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /data/salomonis2/Genomes/Star-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf --limitBAMsortRAM 127779638988
EOF
#for i in *.fastq.gz; do ./STAR.sh $i | bsub; done
 
#for i in *read1.fastq.gz; do ./STAR.sh $i | bsub; done

