#BSUB -W 10:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J altanalyze
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load singularity
singularity run -B /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR/run:/usr/src/app/run \
                --writable /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/altanalyze /usr/src/app/run/bam








