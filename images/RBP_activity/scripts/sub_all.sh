#BSUB -W 24:00
#BSUB -M 196000
#BSUB -n 50
#BSUB -R "span[hosts=1]"
#BSUB -J SRN
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# ./run_stars.py
./run_bbsr.py






































