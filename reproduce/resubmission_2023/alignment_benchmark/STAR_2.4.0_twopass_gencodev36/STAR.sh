#BSUB -W 20:00
#BSUB -M 250G
#BSUB -n 4
#BSUB -R "rusage[mem=250G] span[hosts=1]"
#BSUB -J STAR
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.4.0h
module load parallel

function run_star(){
echo "process $1, first pass"
FASTQ1=/data/salomonis-archive/FASTQs/NCI-R01/TCGA/TCGA-Breast/TCGA-fastq/$1.qsort.read1.fq.gz
FASTQ2=/data/salomonis-archive/FASTQs/NCI-R01/TCGA/TCGA-Breast/TCGA-fastq/$1.qsort.read2.fq.gz
GENOME_DIR=/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/GenomeRef
GENOME=/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/GRCh38.d1.vd1.fa
ROOT_DIR=/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try

STAR --genomeDir ${GENOME_DIR} \
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
     --readFilesCommand gunzip -c \
     --outFileNamePrefix $ROOT_DIR/$1 \
     --outFilterMatchNminOverLread 0.33 \
     --outFilterScoreMinOverLread 0.33 \
     --sjdbOverhang 100 \
     --outSAMstrandField intronMotif \
     --outSAMtype None \
     --outSAMmode None

echo "process $1, build intermediate index for each sample"
mkdir ${GENOME_DIR}_$1
STAR --runMode genomeGenerate \
     --genomeDir ${GENOME_DIR}_$1 \
     --genomeFastaFiles $GENOME \
     --sjdbOverhang 100 \
     --runThreadN 8 \
     --sjdbFileChrStartEnd ${1}SJ.out.tab \
     --outFileNamePrefix $ROOT_DIR/$1

echo "process $1, second pass"
STAR --genomeDir ${GENOME_DIR}_$1 \
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
     --outFileNamePrefix $ROOT_DIR/${1}_second \
     --outFilterMatchNminOverLread 0.33 \
     --outFilterScoreMinOverLread 0.33 \
     --sjdbOverhang 100 \
     --outSAMstrandField intronMotif \
     --outSAMattributes NH HI NM MD AS XS \
     --outSAMunmapped Within \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMheaderHD @HD VN:1.4 


echo "process $1, generate bam index file"
module load samtools
BAM=$ROOT_DIR/${1}_secondAligned.sortedByCoord.out.bam
samtools index $BAM

}



export -f run_star
export TMPDIR=/scratch/ligk2e
cat samples.txt | parallel -P 4 run_star {}
