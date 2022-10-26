#BSUB -W 72:00
#BSUB -M 196G
#BSUB -n 5
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J GTEx_skin
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.4.0h
module load samtools
module load python/2.7.5
module load parallel
module load bedtools

function run(){

    # convert to fastq
    SAMPLE=$1
    ROOT=/data/salomonis-archive/BAMs/NCI-R01/GTEx/newSkin_100922
    BAM=${ROOT}/${SAMPLE}.bam
    samtools sort -n $BAM -o ${ROOT}/Frank/${SAMPLE}.qsort.bam
    bedtools bamtofastq -i ${ROOT}/Frank/${SAMPLE}.qsort.bam -fq ${ROOT}/Frank/${SAMPLE}.1.fastq -fq2 ${ROOT}/Frank/${SAMPLE}.2.fastq
    gzip ${ROOT}/Frank/${SAMPLE}.1.fastq
    gzip ${ROOT}/Frank/${SAMPLE}.2.fastq

    # run STAR
    FASTQ1=${ROOT}/Frank/${SAMPLE}.1.fastq.gz
    FASTQ2=${ROOT}/Frank/${SAMPLE}.2.fastq.gz
    GENOME_DIR=/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/GenomeRef
    GENOME=/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/GRCh38.d1.vd1.fa
    ROOT_DIR=${ROOT}/Frank

    echo "$1 fisrt pass"
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

    echo "$1 intermediate index"
    mkdir ${ROOT_DIR}/GenomeRef_${1}
    STAR --runMode genomeGenerate \
         --genomeDir ${ROOT_DIR}/GenomeRef_${1} \
         --genomeFastaFiles $GENOME \
         --sjdbOverhang 100 \
         --runThreadN 8 \
         --sjdbFileChrStartEnd ${ROOT_DIR}/${1}SJ.out.tab \
         --outFileNamePrefix $ROOT_DIR/$1

    echo "$1 second pass"
    STAR --genomeDir ${ROOT_DIR}/GenomeRef_${1} \
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
         --limitBAMsortRAM 100000000000 \
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

    echo "$1 delete files except bam"
    rm -r ${ROOT_DIR}/GenomeRef_${1}

    rm -r ${ROOT_DIR}/${1}_second_STARtmp
    rm ${ROOT_DIR}/${1}_secondLog.final.out
    rm ${ROOT_DIR}/${1}_secondLog.out
    rm ${ROOT_DIR}/${1}_secondLog.progress.out
    rm ${ROOT_DIR}/${1}_secondSJ.out.tab

    rm ${ROOT_DIR}/${1}Log.final.out
    rm ${ROOT_DIR}/${1}Log.progress.out
    rm ${ROOT_DIR}/${1}SJ.out.tab
    rm ${ROOT_DIR}/${1}Log.out

    rm ${ROOT}/Frank/${SAMPLE}.1.fastq.gz
    rm ${ROOT}/Frank/${SAMPLE}.2.fastq.gz
    rm ${ROOT}/Frank/${SAMPLE}.qsort.bam



    echo "$1 run BAMtoBED"
    BAM=${ROOT_DIR}/${1}_secondAligned.sortedByCoord.out.bam

    python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoJunctionBED.py --i $BAM --species Hs \
           --r /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt

    python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoExonBED.py --i $BAM \
           --r /data/salomonis-archive/BAMs/BMI-ARCHIVES/External-Collaborations/Gladstone/Spindler/STAR-Grch38-results/exonrefdir/Hs.bed --s Hs

    mv ${ROOT_DIR}/${1}_secondAligned.sortedByCoord.out__junction.bed ${ROOT_DIR}/beds
    mv ${ROOT_DIR}/${1}_secondAligned.sortedByCoord.out__intronJunction.bed ${ROOT_DIR}/beds

    echo "$1 remove bams"
    rm ${ROOT_DIR}/${1}_secondAligned.sortedByCoord.out.bam
    rm ${ROOT_DIR}/${1}_secondAligned.sortedByCoord.out.bam.bai
}


export -f run
export TMPDIR=/scratch/ligk2e
cat xaj.txt | parallel -P 5 run {}




