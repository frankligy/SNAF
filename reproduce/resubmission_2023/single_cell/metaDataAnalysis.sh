#BSUB -W 48:00
#BSUB -M 150G
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -J metadataanalysis
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load python/2.7.5

declare -a index_array=(tumor tcell bcell macrophage caf nk endo)
for item in ${index_array[@]}
do 
       python /data/salomonis2/software/AltAnalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no --percentExp 5 \
              --i /data/salomonis-archive/FASTQs/PublicDatasets/scRNA-Seq/aviv_melanoma_Frank/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
              --m /data/salomonis-archive/FASTQs/PublicDatasets/scRNA-Seq/aviv_melanoma_Frank/altanalyze_output/DAS/groups.${item}_other.txt
done

