#BSUB -W 120:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J proliferate
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

module load python/2.7.5

task="original"

### build necessary folder structure
mkdir altanalyze_output
mkdir altanalyze_output/ExpressionInput

### build group file
touch altanalyze_output/ExpressionInput/groups.${task}.txt
cd beds
count=0
for file in *__junction.bed; do 
    stream=$(echo $file | sed 's/__junction.bed/.bed/g')
    if [ $(($count%2)) == 0 ]; then
        stream+='\t1\texp'
    else
        stream+='\t2\tctl'
    fi
    echo -e $stream >> ../altanalyze_output/ExpressionInput/groups.${task}.txt
    ((count+=1))
done
cd ..

### build comp file
touch altanalyze_output/ExpressionInput/comps.${task}.txt
echo -e '1\t2' > altanalyze_output/ExpressionInput/comps.${task}.txt

### run multipath-psi
python /data/salomonis2/software/AltAnalyze-91/AltAnalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
    --bedDir /data/salomonis-archive/FASTQs/PublicDatasets/Bulk-RNASeq/ProliferatingSkinCell/beds \
    --output /data/salomonis-archive/FASTQs/PublicDatasets/Bulk-RNASeq/ProliferatingSkinCell/altanalyze_output \
    --groupdir /data/salomonis-archive/FASTQs/PublicDatasets/Bulk-RNASeq/ProliferatingSkinCell/altanalyze_output/ExpressionInput/groups.${task}.txt \
    --compdir /data/salomonis-archive/FASTQs/PublicDatasets/Bulk-RNASeq/ProliferatingSkinCell/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
    --runGOElite no