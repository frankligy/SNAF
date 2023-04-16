#BSUB -W 60:00
#BSUB -M 196000
#BSUB -n 5
#BSUB -R "span[hosts=1]"
#BSUB -J star
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.6.1


# outSAMstrandFiled, to infer SAM strand from intronMotif.
# sjdbGTFfile, assist alignment by supplying a sjdb.

function run_STAR(){

read file1 file2 <<< $1

STAR --genomeDir /data/salomonis2/Genomes/Star-Index-GRCH38/Grch38-STAR-index \
    --readFilesIn ./fastq/${file1}.fastq.gz ./fastq/${file2}.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ./STAR_output/${file1}_${file2}. \
    --runThreadN 6 \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile /data/salomonis2/Genomes/Star-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf \
    --limitBAMsortRAM 97417671648

return 0
}


# main program
cd ..
module load parallel
export -f run_STAR
export TMPDIR=/scratch/ligk2e    # OR use -tmpdir option for parallel package, default is /tmp
cat xae.txt | parallel -P 5 run_STAR {}



# building the reference

# STAR --runMode genomeGenerate \
#      --genomeDir /data/salomonis2/Genomes/Star-Index-GRCH38/Grch38-STAR-index \
#      --genomeFastaFiles Grch38_r85.all.fa \
#      --runThreadN 8 \
#      --sjdbGTFfile Homo_sapiens.GRCh38.85.gtf \
#      --sjdbOverhang 100

