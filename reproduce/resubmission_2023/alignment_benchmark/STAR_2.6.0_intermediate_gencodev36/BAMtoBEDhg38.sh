#!/bin/bash

BAM=$1
SAMPLE=$(basename $BAM .bam)
DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 2
#BSUB -R "span[ptile=2]"
#BSUB -M 8000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR
module load python/2.7.5
module load samtools

samtools index $BAM

python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoJunctionBED.py --i $BAM --species Hs --r /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt

python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoExonBED.py --i $BAM --r /data/salomonis-archive/BAMs/BMI-ARCHIVES/External-Collaborations/Gladstone/Spindler/STAR-Grch38-results/exonrefdir/Hs.bed --s Hs

EOF
#for i in *.bam; do BAMtoBEDhg38.sh $i | bsub; done
