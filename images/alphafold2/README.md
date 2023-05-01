# How do we run Alphafold2

We first download the community-built singularity image of Alphafold2. Relevant links:

1. https://hub.docker.com/r/catgumag/alphafold
2. https://github.com/dialvarezs/alphafold/blob/main/run_singularity.py 

Then we need to download the data file, we download it to `/database/alphafold/databases`.

With that, the idea now is just to create a singularity command to run, and in the python script, there are some auxiliary function to programatically mount the data, fasta and output directory.

By convention, I create a `fastas` folder to store each fasta file, and I create a `result` folder to generate all predictions.

I want to note this is only for Alphafold2.1.0, it is updating and maybe following the official docker is better.

