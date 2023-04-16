#BSUB -W 24:00
#BSUB -M 196000
#BSUB -n 10
#BSUB -R "span[hosts=1]"
#BSUB -J bam_to_bed
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load python/2.7.5

function run_BAMtoBED() {
 
echo "start to get junction bed"
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoJunctionBED.py --i $1 --species Hs --r /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt
echo "start to get exon bed"
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoExonBED.py --i $1 --r /data/salomonis2/BAMs/BMI-ARCHIVES/External-Collaborations/Gladstone/Spindler/STAR-Grch38-results/exonrefdir/Hs.bed --s Hs  # you can build the reference by youself

return 0
}

#run_BAMtoBED /data/salomonis2/FASTQs/NCI-R01/encode_shrna_kd/K562/STAR_output/ENCFF309THL_ENCFF996JYA.Aligned.sortedByCoord.out.bam

cd ../STAR_output
module load parallel
export -f run_BAMtoBED
export TMPDIR=/scratch/ligk2e
cat ../xaj.txt | parallel -P 10 run_BAMtoBED {}







































