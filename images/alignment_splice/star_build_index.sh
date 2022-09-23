#BSUB -W 10:00
#BSUB -M 250G
#BSUB -n 1
#BSUB -R "rusage[mem=250G] span[hosts=1]"
#BSUB -J STAR
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.6.1

STAR --runMode genomeGenerate \
     --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1 \
     --genomeFastaFiles GRCh38.p13.genome.fa \
     --runThreadN 8 \
     --sjdbGTFfile gencode.v36.annotation.gtf \
     --sjdbOverhang 100







