#BSUB -W 10:00
#BSUB -M 250G
#BSUB -n 4
#BSUB -R "rusage[mem=250G] span[hosts=1]"
#BSUB -J STAR
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


module load STAR/2.6.1
module load parallel

function build_second_index(){

mkdir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1_$1
STAR --runMode genomeGenerate \
     --genomeDir /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/human_gencode_v36_all_scaffold_star_index_v2.6.1_$1 \
     --genomeFastaFiles GRCh38.p13.genome.fa \
     --sjdbOverhang 100 \
     --runThreadN 8 \
     --sjdbFileChrStartEnd ${1}SJ.out.tab \
     --outFileNamePrefix /data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/$1
}

export -f build_second_index
export TMPDIR=/scratch/ligk2e
cat samples.txt | parallel -P 4 build_second_index {}








