#BSUB -W 24:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J GO
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load python/2.7.5

# BioMarkers
mkdir /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/GO_Elite_result_BioMarkers
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/GO_Elite.py --species Hs --mod Ensembl --pval 0.05 --num 3 \
       --input /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/gene_list.txt \
       --output /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/GO_Elite_result_BioMarkers --dataToAnalyze BioMarkers

# GO
mkdir /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/GO_Elite_result_GeneOntology
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/GO_Elite.py --species Hs --mod Ensembl --pval 0.05 --num 3 \
       --input /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/gene_list.txt \
       --output /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/GO_Elite_result_GeneOntology --dataToAnalyze GeneOntology