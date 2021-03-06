#!/bin/bash

cd run
echo "Current folder is "$PWD
echo "BAM file folder is "$1

chmod -R 777 /usr/src/app/altanalyze

# step1: bam to bed


function run_BAMtoBED() {
 
echo "start to get junction bed"
python /usr/src/app/altanalyze/import_scripts/BAMtoJunctionBED.py --i $1 \
    --species Hs --r /usr/src/app/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt

echo "start to get exon bed"
python /usr/src/app/altanalyze/import_scripts/BAMtoExonBED.py --i $1  \
    --r /usr/src/app/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs.bed --s Hs  

return 0
}


cd $1
for file in *.bam; do run_BAMtoBED $file; done
cd ..
mkdir bed
cd $1
for file in *.bed; do mv $file ../bed; done
cd ..

# step2: multipath-psi


# initiate name
task="original"

# build necessary folder structure
mkdir altanalyze_output
mkdir altanalyze_output/ExpressionInput

# build group file
touch altanalyze_output/ExpressionInput/groups.${task}.txt
cd bed
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

# build comp file
touch altanalyze_output/ExpressionInput/comps.${task}.txt
echo -e '1\t2' > altanalyze_output/ExpressionInput/comps.${task}.txt

# run multipath-psi
python /usr/src/app/altanalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
    --bedDir /usr/src/app/run/bed \
    --output /usr/src/app/run/altanalyze_output \
    --groupdir /usr/src/app/run/altanalyze_output/ExpressionInput/groups.${task}.txt \
    --compdir /usr/src/app/run/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
    --runGOElite no

# process the result to get pruned junction count matrix
echo "prune the raw junction count matrix"
python /usr/src/app/prune.py 
