#BSUB -W 20:00
#BSUB -M 196000
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -J alphafold2
#BSUB -q gpu-v100
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


# Dr. kouril
module load singularity/3.7.0
module load python3
unset PYTHONHOME
export CUDA_VISIBLE_DEVICES=1

./run_singularity_mk.py -f ./fastas/novel_DCBLD2.fasta




