#!/bin/bash
#SBATCH --partition=cpu_dev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=4:00:00
#SBATCH --mem=200Gb
#SBATCH --job-name="kallisto"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err

module load parallel
module load kallisto/0.44.0

function run() {
    full=$(echo $1 | sed 's/-/\//' | sed 's/_secondAligned.sortedByCoord.out//')
    sample=$(echo $1 | sed 's/_secondAligned.sortedByCoord.out//')
    mkdir process/${sample}
    kallisto quant -i transcripts.idx -o process/${sample} -b 100 --single -l 180 -s 20 /gpfs/data/yarmarkovichlab/senescence/RNA_seq_fastq/${full}_R1_001.fastq.gz
    ./t2g_tpm.py --input process/${sample}/abundance.tsv

}

export -f run
export TMPDIR=/gpfs/scratch/lig08
cat ../samples.txt | parallel -P 30 run {}