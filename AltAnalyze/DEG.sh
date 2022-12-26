#!/bin/bash

python /usr/src/AltAnalyze/stats_scripts/metaDataAnalysis.py --p RNASeq --s Hs --adjp yes --pval 1 --f 1 \
       --i $1 \
       --m /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/groups.txt 