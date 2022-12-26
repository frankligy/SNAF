#!/bin/bash


function fastq_to_junction(){
     # build index
     ./STAR_index.sh

     # pass the mode and node
     ./STAR_align.sh $1 $2

     # multipath-psi
     ./AltAnalyze.sh 
}

function DEG(){
     ./DEG.sh $1
}

function GO_Elite(){
     ./GO_Elite.sh $1
}


while getopts "k:mni" flag; do
     case $flag in
     k)
     kind=$OPTARG
     ;;
     m)
     mode=$OPTARG
     ;;
     n)
     node=$OPTARG
     ;;
     i)
     input=$OPTARG
     ;;
     \?)
     echo "misformated arguments"
     exit 1
     ;;
     esac
done

if [ kind == 'count' ]
     then fastq_to_junction ${mode} ${node}
elif [ kind == 'DEG' ]
     then DEG ${input}
elif [ kind == 'GO_Elite' ]
     then GO_Elite ${input}