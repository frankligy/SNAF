#!/bin/bash

FASTQ1=$1
FASTQ2=${FASTQ1/read1/read2}
SAMPLE=$(basename $FASTQ1 .fq.gz)

DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 8
#BSUB -R "span[ptile=4]"
#BSUB -M 128000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR
module load STAR/2.7.4

STAR --runMode genomeGenerate --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR2pass-GDCRef/GenomeRef/SJV36 --genomeFastaFiles /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR2pass-GDCRef/GenomeRef/GRCh38.d1.vd1.fa --sjdbGTFfile /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR2pass-GDCRef/GenomeRef/gencode.v36.annotation.gtf --runThreadN 8 --sjdbOverhang 100 --sjdbFileChrStartEnd /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR2pass/SJ_out/*SJ.out.tab

EOF
#for i in *.fastq.gz; do ./STAR.sh $i | bsub; done
 
#for i in *read1.fastq.gz; do ./STAR.sh $i | bsub; done

