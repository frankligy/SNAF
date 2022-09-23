# Testing the proper ways to run STAR for junction detection

### Resources

1. [STAR manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
2. [GDC TCGA manual](https://docs.gdc.cancer.gov/Data/PDF/Data_UG.pdf)
3. [SAM manual](https://en.wikipedia.org/wiki/SAM_(file_format))

### Four testing fastq files converted from TCGA BAM

```
Root path: 
/data/salomonis-archive/FASTQs/NCI-R01/TCGA/TCGA-Breast/TCGA-fastq

[1] sample1 (read1 and read2)
TCGA-A1-A0SI-01A-11R-A144-07.qsort.read1.fq.gz

[2] sample2 (read1 and read2)
TCGA-A1-A0SM-01A-11R-A084-07.qsort.read1.fq.gz

[3] sample3 (read1 and read2)
TCGA-A1-A0SO-01A-22R-A084-07.qsort.read1.fq.gz

[4] sample4 (read1 and read2)
TCGA-A1-A0SN-01A-11R-A144-07.qsort.read1.fq.gz
```

### TCGA BAM (ground-state truth)

```
Root path:
/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BAM

[1] sample1
TCGA-A1-A0SI-01A-11R-A144-07.bam

[2] sample2
TCGA-A1-A0SM-01A-11R-A084-07.bam

[3] sample3
TCGA-A1-A0SO-01A-22R-A084-07.bam

[4] sample4
TCGA-A1-A0SN-01A-11R-A144-07.bam
```

### Pipeline

```bash
# build initial index
module load STAR/2.6.1
STAR --runMode genomeGenerate \
     --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1 \
     --genomeFastaFiles GRCh38.p13.genome.fa \
     --runThreadN 8 \
     --sjdbGTFfile gencode.v36.annotation.gtf \
     --sjdbOverhang 100
```

**Notes**

1. Hard to figure out STAR version they used
2. using gencode v36 as they mentioned
3. include other scaffold as suggested in STAR manual


```bash
# pass1
STAR --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1 \
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
     --outFileNamePrefix $DIR/$SAMPLE \
     --outFilterMatchNminOverLread 0.33
     --outFilterScoreMinOverLread 0.33
     --sjdbOverhang 100
     --outSAMstrandField intronMotif
     --outSAMtype None
     --outSAMmode None
```

**Notes**

1. The only difference is add `--outFileNamePrefix`

```bash
# build second index
STAR --runMode genomeGenerate \
     --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1_$1 \
     --genomeFastaFiles GRCh38.p13.genome.fa \
     --sjdbOverhang 100 \
     --runThreadN 8 \
     --sjdbFileChrStartEnd ${1}SJ.out.tab \
     --outFileNamePrefix /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/$1
```

**Notes**

1. Add `outFileNamePrefix` so that each sample will generate its own prefix_STARtmp folder


```bash
# pass2
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
     --outFileNamePrefix $DIR/$SAMPLE \
     --outFilterMatchNminOverLread 0.33 \
     --outFilterScoreMinOverLread 0.33 \
     --sjdbOverhang 100 \
     --outSAMstrandField intronMotif \
     --outSAMattributes NH HI NM MD AS XS \
     --outSAMunmapped Within \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMheaderHD @HD VN:1.4 
```

**Notes**

1. Still need to run samtools to generate bam.bai

### Testing 

```
redacted
```
