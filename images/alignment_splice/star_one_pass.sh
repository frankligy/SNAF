#BSUB -W 10:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J STAR_splice
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# align to the pre-built index
module load STAR/2.6.1
FASTQ1=/data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/SRR12086659_1.fastq
FASTQ2=/data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/SRR12086659_2.fastq
STAR --genomeDir /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/mouse_genome_index \
     --readFilesIn $FASTQ1 $FASTQ2 \
     --outFileNamePrefix SRR12086659 \
     --runThreadN 6 \
     --outSAMstrandField intronMotif \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbGTFfile /data/salomonis2/LabFiles/Frank-Li/upstream/STAR_gene_count/mouse/gencode.vM30.chr_patch_hapl_scaff.annotation.gtf \
     --limitBAMsortRAM 97417671648








