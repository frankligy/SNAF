#BSUB -W 10:00
#BSUB -M 196G
#BSUB -n 1
#BSUB -R "rusage[mem=196G] span[hosts=1]"
#BSUB -J TCGA_melanoma
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# run snaf
module load singularity
./analysis_new.py

# # run maxquant
# function run(){
#     cd /data/salomonis-archive/MS/melanoma/raw/${1}
#     mono /data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe /data/salomonis-archive/MS/melanoma/raw/${1}/mqpar.xml
# }

# module load mono/5.20.1
# module load dotnet
# export PATH=$PATH:​/usr/local/mono/5.20.1/bin

# while read line
#     do run $line
# done < ms_patients.txt



