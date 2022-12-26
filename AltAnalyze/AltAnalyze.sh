#!/bin/bash


# step1: bam to bed
function run_BAMtoBED() {
 
echo "start to get junction bed"
python /usr/src/altanalyze/import_scripts/BAMtoJunctionBED.py --i $1 \
    --species Hs --r /usr/src/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt

echo "start to get exon bed"
python /usr/src/altanalyze/import_scripts/BAMtoExonBED.py --i $1  \
    --r /usr/src/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs.bed --s Hs  

return 0
}

### collect for bam file name for parallelization
cd /mnt/bams
for file in *.bam; do echo $file; done > ../samples.txt 

### start to run
export -f run_BAMtoBED
export TMPDIR=/tmp
cat ../samples.txt | parallel -P $2 run_BAMtoBED {}

### move bed files to bed folder
cd ..
mkdir beds
cd $1
for file in *.bed; do mv $file ../beds; done
cd ..

# step2: multipath-psi
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
python /usr/src/altanalyze/AltAnalyze.py --species Hs --platform RNASeq --version EnsMart91 \
    --bedDir /mnt/beds \
    --output /mnt/altanalyze_output \
    --groupdir /mnt/altanalyze_output/ExpressionInput/groups.${task}.txt \
    --compdir /mnt/altanalyze_output/ExpressionInput/comps.${task}.txt --expname ${task} \
    --runGOElite no

# step3: process count matrix to only contain PSI junctions
echo "prune the raw junction count matrix"
python /usr/src/prune.py 
