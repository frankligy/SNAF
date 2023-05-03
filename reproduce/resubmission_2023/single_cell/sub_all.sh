#BSUB -W 10:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J aviv
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


# ./analysis.py

# concat fastq
cd ./typing_fastq
declare -a index_array=(53 58 60 67 71 72 74 75 78 79 80 81 82 84 88 89 94)
for sample in ${index_array[@]}
do 
    cd ./${sample}
    for file in *.fastq.gz; do gunzip $file; done
    cat *R1_001.fastq > ../${sample}_combine_R1.fastq
    cat *R2_001.fastq > ../${sample}_combine_R2.fastq
    cd ..
done







