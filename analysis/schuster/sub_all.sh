#BSUB -W 20:00
#BSUB -M 100000
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# download
# module load sratoolkit/2.10.4
# fasterq-dump -e 20 SRR5933736
# fasterq-dump -e 20 SRR5933737

# run STAR
# module load STAR/2.6.1
# declare -a index_array
# index_array+=(SRR5933736 SRR5933737)
# for srr in ${index_array[@]}
# do
#     STAR --genomeDir /data/salomonis2/Genomes/Star-Index-GRCH38/Grch38-STAR-index \
#         --readFilesIn $(pwd)/$srr.fastq \
#         --outFileNamePrefix $(pwd)/$srr. \
#         --runThreadN 4 \
#         --outSAMstrandField intronMotif \
#         --outSAMtype BAM SortedByCoordinate \
#         --sjdbGTFfile /data/salomonis2/Genomes/Star-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf \
#         --limitBAMsortRAM 97417671648
# done

# run bam_to_bed
# module load python/2.7.5
# declare -a index_array
# index_array+=(SRR5933736.Aligned.sortedByCoord.out.bam SRR5933737.Aligned.sortedByCoord.out.bam)
# for bam in ${index_array[@]}
# do 
#     echo "start to get junction bed"
#     python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoJunctionBED.py --i $bam --species Hs --r /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt
#     echo "start to get exon bed"
#     python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/import_scripts/BAMtoExonBED.py --i $bam --r /data/salomonis2/BAMs/BMI-ARCHIVES/External-Collaborations/Gladstone/Spindler/STAR-Grch38-results/exonrefdir/Hs.bed --s Hs  
# done


# run altanalyze with existing bed files
module load samtools 
module load python/2.7.5
module load R

task="original"
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
    --bedDir ./RNA/junction --output ./RNA/altanalyze_output --groupdir ./altanalyze_output/ExpressionInput/groups.${task}.txt \
    --compdir ./altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
    --runGOElite no


