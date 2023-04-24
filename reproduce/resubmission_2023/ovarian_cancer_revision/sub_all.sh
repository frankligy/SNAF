#BSUB -W 12:00
#BSUB -M 100G
#BSUB -n 20
#BSUB -R "rusage[mem=100G] span[hosts=1]"
#BSUB -J ovarian
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# # run SNAF
# ./analysis.py

# run maxQuant
SAMPLE=OvCa48
module load mono/5.20.1
module load dotnet
export PATH=$PATH:​/usr/local/mono/5.20.1/bin
mono /data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe /data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ovarian/MS/${SAMPLE}/mqpar.xml

