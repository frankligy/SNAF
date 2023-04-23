#BSUB -W 10:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J altanalyze
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load singularity
cd bam
singularity run -B $PWD:/usr/src/app/run --writable ../altanalyze /usr/src/app/run/bam








