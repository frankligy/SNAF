#BSUB -W 10:00
#BSUB -M 250G
#BSUB -n 4
#BSUB -R "rusage[mem=250G] span[hosts=1]"
#BSUB -J STAR
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.6.1
module load parallel

function run_star(){
echo "process $1"
FASTQ1=/data/salomonis-archive/FASTQs/NCI-R01/TCGA/TCGA-Breast/TCGA-fastq/$1.qsort.read1.fq.gz
FASTQ2=/data/salomonis-archive/FASTQs/NCI-R01/TCGA/TCGA-Breast/TCGA-fastq/$1.qsort.read2.fq.gz
SAMPLE=$1
STAR --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1_$1 \
     --readFilesIn $FASTQ1 $FASTQ2 \
     --runThreadN 8 \
     --outFilterMultimapScoreRange 1 \
     --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 10 \
     --alignIntronMax 500000 \
     --alignMatesGapMax 1000000 \
     --sjdbScore 2 \
     --alignSJDBoverhangMin 1 \
     --genomeLoad NoSharedMemory \
     --limitBAMsortRAM 0 \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/${1}_second \
     --outFilterMatchNminOverLread 0.33 \
     --outFilterScoreMinOverLread 0.33 \
     --sjdbOverhang 100 \
     --outSAMstrandField intronMotif \
     --outSAMattributes NH HI NM MD AS XS \
     --outSAMunmapped Within \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMheaderHD @HD VN:1.4 
}

export -f run_star
export TMPDIR=/scratch/ligk2e
cat samples.txt | parallel -P 4 run_star {}









