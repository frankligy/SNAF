#BSUB -W 48:00
#BSUB -M 150000
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -J metadataanalysis
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load python/2.7.5

function program_shRNA_K562 (){
# copy the group file and rename it
bn=$(basename $1)  # groups.YBX3_shRNA_K562.txt
arrbn=(${bn//./ })
id=${arrbn[1]}   # YBX3_shRNA_K562
arrid=(${id//_/ })
sf=${arrid[0]} # YBX3
cp $1 /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/${bn} /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# running metaDataAnaysis
python /data/salomonis2/software/AltAnalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
       --i /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
       --m /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# copy the psi file to final location and rename it
cp /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562/${sf}.txt
# remove the intermediate results
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/ExpressionProfiles
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/PValues
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.pdf
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.png
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/top50
}

function program_shRNA_HepG2 (){
# copy the group file and rename it
bn=$(basename $1)  # groups.YBX3_shRNA_K562.txt
arrbn=(${bn//./ })
id=${arrbn[1]}   # YBX3_shRNA_K562
arrid=(${id//_/ })
sf=${arrid[0]} # YBX3
cp $1 /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/${bn} /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# running metaDataAnaysis
python /data/salomonis2/software/AltAnalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
       --i /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
       --m /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# copy the psi file to final location and rename it
cp /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2/${sf}.txt
# remove the intermediate results
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/ExpressionProfiles
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/PValues
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.pdf
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.png
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/top50
}

function program_CRISPR_K562 (){
# copy the group file and rename it
bn=$(basename $1)  # groups.YBX3_shRNA_K562.txt
arrbn=(${bn//./ })
id=${arrbn[1]}   # YBX3_shRNA_K562
arrid=(${id//_/ })
sf=${arrid[0]} # YBX3
cp $1 /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/${bn} /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# running metaDataAnaysis
python /data/salomonis2/software/AltAnalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
       --i /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
       --m /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# copy the psi file to final location and rename it
cp /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_K562
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_K562/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_K562/${sf}.txt
# remove the intermediate results
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/ExpressionProfiles
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/PValues
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.pdf
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.png
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/top50
}

function program_CRISPR_HepG2 (){
# copy the group file and rename it
bn=$(basename $1)  # groups.YBX3_shRNA_K562.txt
arrbn=(${bn//./ })
id=${arrbn[1]}   # YBX3_shRNA_K562
arrid=(${id//_/ })
sf=${arrid[0]} # YBX3
cp $1 /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/${bn} /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# running metaDataAnaysis
python /data/salomonis2/software/AltAnalyze/stats_scripts/metaDataAnalysis.py --p PSI --dPSI 0 --pval 1 --adjp no \
       --i /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt \
       --m /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
# copy the psi file to final location and rename it
cp /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_HepG2
mv /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_HepG2/PSI.exp_vs_control.txt \
   /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_HepG2/${sf}.txt
# remove the intermediate results
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Events-dPSI_0.0_rawp
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/ExpressionProfiles
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/PValues
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.pdf
rm /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/*.png
rm -r /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/top50
}


mkdir /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562
for file in /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_combine/groups.*_shRNA_K562.txt; do
program_shRNA_K562 $file
done

mkdir /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2
for file in /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_combine/groups.*_shRNA_HepG2.txt; do
program_shRNA_HepG2 $file
done

mkdir /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_K562
for file in /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_combine/groups.*_CRISPR_K562.txt; do
program_CRISPR_K562 $file
done

mkdir /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/CRISPR_HepG2
for file in /data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_combine/groups.*_CRISPR_HepG2.txt; do
program_CRISPR_HepG2 $file
done












