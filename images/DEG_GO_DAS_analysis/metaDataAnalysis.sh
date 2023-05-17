#BSUB -W 24:00
#BSUB -M 150G
#BSUB -n 1
#BSUB -R "rusage[mem=150G] span[hosts=1]"
#BSUB -J metadataanalysis
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load python/2.7.5


# running metaDataAnaysis
python /data/salomonis2/software/AltAnalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
       --i /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
       --m /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/survival/groups.txt


