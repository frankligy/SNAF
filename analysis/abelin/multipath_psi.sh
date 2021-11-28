#BSUB -W 10:00
#BSUB -M 196000
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -J multipath_psi
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load samtools 
module load python/2.7.5
module load R

task="original"
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
    --bedDir ./junction --output ./altanalyze_output --groupdir ./altanalyze_output/ExpressionInput/groups.${task}.txt \
    --compdir ./altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
    --runGOElite no