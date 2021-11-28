#BSUB -W 10:00
#BSUB -M 196000
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -J star_align
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.6.1


# outSAMstrandFiled, to infer SAM strand from intronMotif.
# sjdbGTFfile, assist alignment by supplying a sjdb.

function run_STAR(){

read file1 file2 <<< $1

STAR --genomeDir /data/salomonis2/Genomes/Star-Index-GRCH38/Grch38-STAR-index \
    --readFilesIn ./${file1}.fastq.gz ./${file2}.fastq.gz \
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
module load parallel
export -f run_STAR
cat pair_info.txt | parallel -P 4 run_STAR {}