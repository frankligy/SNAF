# How do we run Alphafold2

We first download the community-built singularity image of Alphafold2. Relevant links:

1. https://hub.docker.com/r/catgumag/alphafold
2. https://github.com/dialvarezs/alphafold/blob/main/run_singularity.py 

Then we need to download the data file, we download it to `/database/alphafold/databases`.

With that, the idea now is just to create a singularity command to run, and in the python script, there are some auxiliary function to programatically mount the data, fasta and output directory.

By convention, I create a `fastas` folder to store each fasta file, and I create a `result` folder to generate all predictions.

I want to note this is only for Alphafold2.1.0, it is updating and maybe following the official docker is better.

Below is an example running command:

```bash
singularity exec --nv --bind /data/salomonis2/LabFiles/Frank-Li/neoantigen/alphafold/fastas:/mnt/fasta_path_0:ro,/database/alphafold/databases/uniref90:/mnt/uniref90_database_path:ro,/database/alphafold/databases/mgnify:/mnt/mgnify_database_path:ro,/database/alphafold:/mnt/data_dir:ro,/database/alphafold/databases/pdb_mmcif:/mnt/template_mmcif_dir:ro,/database/alphafold/databases/pdb_mmcif:/mnt/obsolete_pdbs_path:ro,/database/alphafold/databases/pdb70:/mnt/pdb70_database_path:ro,/database/alphafold/databases/uniclust30/uniclust30_2018_08:/mnt/uniclust30_database_path:ro,/database/alphafold/databases/bfd:/mnt/bfd_database_path:ro,/data/salomonis2/LabFiles/Frank-Li/neoantigen/alphafold:/mnt/output:rw --env="NVIDIA_VISIBLE_DEVICES=all" --env="TF_FORCE_UNIFIED_MEMORY=1" --env="XLA_PYTHON_CLIENT_MEM_FRACTION=4.0" --env="MAX_CPUS=8" /data/salomonis2/LabFiles/Frank-Li/neoantigen/alphafold/alphafold.sif /app/run_alphafold.sh --fasta_paths=/mnt/fasta_path_0/ref_IGSF11.fasta --uniref90_database_path=/mnt/uniref90_database_path/uniref90.fasta --mgnify_database_path=/mnt/mgnify_database_path/mgy_clusters_2018_12.fa --data_dir=/mnt/data_dir/databases --template_mmcif_dir=/mnt/template_mmcif_dir/mmcif_files --obsolete_pdbs_path=/mnt/obsolete_pdbs_path/obsolete.dat --pdb70_database_path=/mnt/pdb70_database_path/pdb70 --uniclust30_database_path=/mnt/uniclust30_database_path/uniclust30_2018_08 --bfd_database_path=/mnt/bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --output_dir=/mnt/output/results --max_template_date=2023-05-01 --db_preset=full_dbs --model_preset=monomer --benchmark=False --use_precomputed_msas=False --logtostderr
```

