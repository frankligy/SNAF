#!/bin/bash



function pair_end_run(){

    # pass SRR3333 to this function, require your file to be sample.R1.fastq.gz and sample.R2.fastq.gz
    SAMPLE=$1
    FASTQ_DIR=/mnt/fastqs
    FASTQ1=${FASTQ_DIR}/${SAMPLE}.R1.fastq.gz
    FASTQ1=${FASTQ_DIR}/${SAMPLE}.R2.fastq.gz
    GENOME_DIR=/tmp/GenomeRef
    GENOME=/tmp/GRCh38.d1.vd1.fa
    ROOT_DIR=/mnt
    mkdir /mnt/bams

    echo "${SAMPLE} fisrt pass"
    STAR --genomeDir ${GENOME_DIR} \
         --readFilesIn ${FASTQ1} ${FASTQ2} \
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
         --outFileNamePrefix ${ROOT_DIR}/${SAMPLE} \
         --outFilterMatchNminOverLread 0.33 \
         --outFilterScoreMinOverLread 0.33 \
         --sjdbOverhang 100 \
         --outSAMstrandField intronMotif \
         --outSAMtype None \
         --outSAMmode None

    echo "${SAMPLE} intermediate index"
    mkdir ${ROOT_DIR}/GenomeRef_${SAMPLE}
    STAR --runMode genomeGenerate \
         --genomeDir ${ROOT_DIR}/GenomeRef_${SAMPLE} \
         --genomeFastaFiles ${GENOME} \
         --sjdbOverhang 100 \
         --runThreadN 8 \
         --sjdbFileChrStartEnd ${ROOT_DIR}/${SAMPLE}SJ.out.tab \
         --outFileNamePrefix ${ROOT_DIR}/${SAMPLE}

    echo "${SAMPLE} second pass"
    STAR --genomeDir ${ROOT_DIR}/GenomeRef_${1} \
         --readFilesIn ${FASTQ1} ${FASTQ2} \
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
         --outFileNamePrefix ${ROOT_DIR}/${SAMPLE}_second \
         --outFilterMatchNminOverLread 0.33 \
         --outFilterScoreMinOverLread 0.33 \
         --sjdbOverhang 100 \
         --outSAMstrandField intronMotif \
         --outSAMattributes NH HI NM MD AS XS \
         --outSAMunmapped Within \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMheaderHD @HD VN:1.4 

    echo "${SAMPLE} delete files except bam and move bam to bam folder"
    rm -r ${ROOT_DIR}/GenomeRef_${SAMPLE}

    rm -r ${ROOT_DIR}/${SAMPLE}_second_STARtmp
    rm ${ROOT_DIR}/${SAMPLE}_secondLog.final.out
    rm ${ROOT_DIR}/${SAMPLE}_secondLog.out
    rm ${ROOT_DIR}/${SAMPLE}_secondLog.progress.out
    rm ${ROOT_DIR}/${SAMPLE}_secondSJ.out.tab

    rm ${ROOT_DIR}/${SAMPLE}Log.final.out
    rm ${ROOT_DIR}/${SAMPLE}Log.progress.out
    rm ${ROOT_DIR}/${SAMPLE}SJ.out.tab
    rm ${ROOT_DIR}/${SAMPLE}Log.out

    mv ${ROOT_DIR}/${SAMPLE}_secondAligned.sortedByCoord.out.bam /mnt/bams
}

function single_end_run(){

    # pass SRR3333 to this function, require your file to be sample.fastq.gz
    SAMPLE=$1
    FASTQ_DIR=/mnt/fastqs
    FASTQ=${FASTQ_DIR}/${SAMPLE}.fastq.gz
    GENOME_DIR=/tmp/GenomeRef
    GENOME=/tmp/GRCh38.d1.vd1.fa
    ROOT_DIR=/mnt
    mkdir /mnt/bams

    echo "${SAMPLE} fisrt pass"
    STAR --genomeDir ${GENOME_DIR} \
         --readFilesIn ${FASTQ} \
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
         --outFileNamePrefix ${ROOT_DIR}/${SAMPLE} \
         --outFilterMatchNminOverLread 0.33 \
         --outFilterScoreMinOverLread 0.33 \
         --sjdbOverhang 100 \
         --outSAMstrandField intronMotif \
         --outSAMtype None \
         --outSAMmode None

    echo "${SAMPLE} intermediate index"
    mkdir ${ROOT_DIR}/GenomeRef_${SAMPLE}
    STAR --runMode genomeGenerate \
         --genomeDir ${ROOT_DIR}/GenomeRef_${SAMPLE} \
         --genomeFastaFiles ${GENOME} \
         --sjdbOverhang 100 \
         --runThreadN 8 \
         --sjdbFileChrStartEnd ${ROOT_DIR}/${SAMPLE}SJ.out.tab \
         --outFileNamePrefix ${ROOT_DIR}/${SAMPLE}

    echo "${SAMPLE} second pass"
    STAR --genomeDir ${ROOT_DIR}/GenomeRef_${1} \
         --readFilesIn ${FASTQ} \
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
         --outFileNamePrefix ${ROOT_DIR}/${SAMPLE}_second \
         --outFilterMatchNminOverLread 0.33 \
         --outFilterScoreMinOverLread 0.33 \
         --sjdbOverhang 100 \
         --outSAMstrandField intronMotif \
         --outSAMattributes NH HI NM MD AS XS \
         --outSAMunmapped Within \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMheaderHD @HD VN:1.4 

    echo "${SAMPLE} delete files except bam and move bam to bam folder"
    rm -r ${ROOT_DIR}/GenomeRef_${SAMPLE}

    rm -r ${ROOT_DIR}/${SAMPLE}_second_STARtmp
    rm ${ROOT_DIR}/${SAMPLE}_secondLog.final.out
    rm ${ROOT_DIR}/${SAMPLE}_secondLog.out
    rm ${ROOT_DIR}/${SAMPLE}_secondLog.progress.out
    rm ${ROOT_DIR}/${SAMPLE}_secondSJ.out.tab

    rm ${ROOT_DIR}/${SAMPLE}Log.final.out
    rm ${ROOT_DIR}/${SAMPLE}Log.progress.out
    rm ${ROOT_DIR}/${SAMPLE}SJ.out.tab
    rm ${ROOT_DIR}/${SAMPLE}Log.out

    mv ${ROOT_DIR}/${SAMPLE}_secondAligned.sortedByCoord.out.bam /mnt/bams
}


mode=$1
node=$2
export -f single_end_run
export -f pair_end_run
export TMPDIR=/tmp

if [ ${mode} == "single" ]
then
    cd /mnt/fastqs
    for file in *.fastq.gz; do echo $(basename -s .fastq.gz $file); done > /mnt/samples.txt 
    cd /mnt
    cat samples.txt | parallel -P ${node} single_end_run {}
elif [ ${mode} == "pair" ]
then
    cd /mnt/fastqs
    for file in *.R1.fastq.gz; do echo $(basename -s .R1.fastq.gz $file); done > /mnt/samples.txt
    cd /mnt
    cat samples.txt | parallel -P ${node} pair_end_run {}
else
    echo "mode either be single or pair"
fi



