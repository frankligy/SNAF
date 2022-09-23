#BSUB -W 10:00
#BSUB -M 250G
#BSUB -n 4
#BSUB -R "rusage[mem=250G] span[hosts=1]"
#BSUB -J STAR
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load samtools
module load parallel

function generate_bai(){
echo "process $1"
BAM=/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/${1}_secondAligned.sortedByCoord.out.bam
samtools index $BAM
}

export -f generate_bai
export TMPDIR=/scratch/ligk2e
cat samples.txt | parallel -P 4 generate_bai {}









