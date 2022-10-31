#BSUB -W 24:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J DEG
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load python/2.7.5

python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/stats_scripts/metaDataAnalysis.py --p RNASeq --s Hs --adjp yes --pval 0.05 --f 2 \
       --i /data/salomonis2/NCI-R01/TCGA-SKCM/TCGA_SKCM_hg38/ExpressionInput/exp.TCGA-SKCM.txt \
       --m /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/groups.txt 

# the result will be generated in the ExpressionInput folder