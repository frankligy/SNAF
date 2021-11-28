#BSUB -W 10:00
#BSUB -M 32000
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -J sratool
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err



# # download
# module load sratoolkit/2.10.4
# cut -f 1 srr_list.txt | xargs -L 1 -P 4 fasterq-dump 

# compress
for file in *.fastq; do gzip $file; done

