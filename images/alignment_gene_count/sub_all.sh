#BSUB -W 10:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J STAR_mouse_gene_count
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# # download
# module load sratoolkit/2.10.4
# fasterq-dump -e 20 SRR12086659

# build index
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/GRCm39.genome.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.chr_patch_hapl_scaff.annotation.gtf.gz
# gunzip -c GRCm39.genome.fa.gz > GRCm39.genome.fa
# gunzip -c gencode.vM30.chr_patch_hapl_scaff.annotation.gtf.gz > gencode.vM30.chr_patch_hapl_scaff.annotation.gtf
# module load STAR/2.6.1
# mkdir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index
# STAR --runMode genomeGenerate \
#      --genomeDir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index \
#      --genomeFastaFiles GRCm39.genome.fa \
#      --runThreadN 8 \
#      --sjdbGTFfile gencode.vM30.chr_patch_hapl_scaff.annotation.gtf \
#      --sjdbOverhang 40

# align
module load STAR/2.6.1
STAR --runThreadN 8 \
     --genomeDir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index \
     --readFilesIn test.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFileNamePrefix testing_









