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

STAR --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR2pass-GDCRef/GenomeRef/SJV36 --readFilesIn $FASTQ1 $FASTQ2 --readFilesCommand gunzip -c --outFileNamePrefix $DIR/$SAMPLE --runThreadN 4 --outSAMunmapped Within --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMattributes NH HI NM MD AS XS



EOF
#for i in *.fastq.gz; do ./STAR.sh $i | bsub; done
 
#for i in *read1.fq.gz; do ./STAR2pass_step3.sh $i | bsub; done

