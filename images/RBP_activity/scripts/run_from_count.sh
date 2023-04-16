#BSUB -W 24:00
#BSUB -M 100000
#BSUB -n 1
#BSUB -J ExpressionBuilder
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


# running ExpressionBuilder.py
module load AltAnalyze

mkdir /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bc
cd /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bc
cp -R /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/AltDatabase ./
mkdir ExpressionInput
cd ./ExpressionInput
cp /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.limma.processed.txt ./
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/ExpressionBuilder.py --species Hs --platform RNASeq \
    --i "/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bc" --a unbiased \
    --additional "/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bc/ExpressionInput/counts.original.limma.processed.txt"

# Continue to run AugmentEventAnnotations.py
module unload AltAnalyze
module load samtools 
module load python/2.7.5
module load R
module load Bioconductor/3.4-R_3.3.1

cd /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bcAltResults  
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/AugmentEventAnnotations.py --i /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bcAltResults/AlternativeOutput --species Hs

