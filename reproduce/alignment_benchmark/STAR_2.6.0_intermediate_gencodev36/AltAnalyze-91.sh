#!/bin/bash

BAM=$1
SAMPLE=$(basename $BAM .txt)
DIR=$(pwd)

cat <<EOF
#BSUB -L /bin/bash
#BSUB -W 60:00
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -M 96000
#BSUB -e $DIR/logs/%J.err
#BSUB -o $DIR/logs/%J.out
#BSUB -J $SAMPLE

cd $DIR

mkdir -p logs

module load samtools 
module load python/2.7.5
module load R


python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 --bedDir $DIR --groupdir $DIR/ExpressionInput/groups.test.txt --compdir $DIR/ExpressionInput/comps.test.txt --expname test --runGOElite no


EOF
#./AltAnalyze-91.sh | bsub