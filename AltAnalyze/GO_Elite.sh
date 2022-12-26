#!/bin/bash

# BioMarkers
mkdir /mnt/GO_Elite_result_BioMarkers
python /usr/src/AltAnalyze/GO_Elite.py --species Hs --mod Ensembl --pval 0.05 --num 3 \
       --input $1 \
       --output /mnt/GO_Elite_result_BioMarkers --dataToAnalyze BioMarkers

# GO
mkdir /mnt/GO_Elite_result_GeneOntology
python /usr/src/AltAnalyze/GO_Elite.py --species Hs --mod Ensembl --pval 0.05 --num 3 \
       --input $1 \
       --output /mnt/GO_Elite_result_GeneOntology --dataToAnalyze GeneOntology