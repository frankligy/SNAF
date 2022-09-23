## 1. Download a pair-end mouse FASTQ file from GEO

```bash
module load sratoolkit/2.10.4
fasterq-dump -e 20 SRR12086659
```

## 2. Align using STAR (build Mouse Index)

```bash
# genome sequence
wget ftps://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.genome.fa.gz

# genome annotation
wget ftps://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.chr_patch_hapl_scaff.annotation.gtf.gz

# decompress
gunzip -c GRCm39.genome.fa.gz > GRCm39.genome.fa
gunzip -c gencode.vM30.chr_patch_hapl_scaff.annotation.gtf.gz > gencode.vM30.chr_patch_hapl_scaff.annotation.gtf

# build STAR index
module load STAR/2.6.1
mkdir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index
STAR --runMode genomeGenerate \
     --genomeDir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index \
     --genomeFastaFiles GRCm39.genome.fa \
     --runThreadN 8 \
     --sjdbGTFfile gencode.vM30.chr_patch_hapl_scaff.annotation.gtf \
     --sjdbOverhang 40
```

## 3. Align with STAR (align and count gene read)

```bash
module load STAR/2.6.1
STAR --runThreadN 8 \
     --genomeDir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index \
     --readFilesIn SRR12086659_1.fastq SRR12086659_2.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFileNamePrefix testing_
```