#BSUB -W 20:00
#BSUB -M 100000
#BSUB -n 1
#BSUB -J download
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


# cd ./bed
# cat ../bed_url_double_quoted.txt | xargs -L 1 curl -O -J -L 

module load anaconda3
source activate ./idr_new_env
cd ./bed
../assist.py





