#BSUB -W 20:00
#BSUB -M 196000
#BSUB -n 6
#BSUB -J OptiType
#BSUB -R "span[hosts=1]"
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load singularity/3.1.0

# the sif is obtained by singularity build my_software.sif docker://fred2/optitype

function program() {

echo $1

sample=$1
mkdir $(pwd)/process/${sample}

echo 'copying files, gunzip and renaming'
read1=/data/salomonis-archive/FASTQs/NCI-R01/Wise-Draper/${sample}.R1.fastq.gz
read2=/data/salomonis-archive/FASTQs/NCI-R01/Wise-Draper/${sample}.R2.fastq.gz
cd $(pwd)/process/${sample}
cp ${read1} ./
cp ${read2} ./
gunzip ${sample}.R1.fastq.gz
gunzip ${sample}.R2.fastq.gz

echo 'running Optitype through singularity'
singularity run -B $(pwd):/mnt ../../my_software.sif -i ${sample}.R1.fastq ${sample}.R2.fastq --rna -v -o /mnt

echo 'remove fastq files'
rm ${sample}.R1.fastq
rm ${sample}.R2.fastq

}

module load parallel
export -f program
export TMPDIR=/scratch/ligk2e
cat xag.txt | parallel -P 6 program {}












