Installation
===============

SNAF consists of two major analysis components. The first is unbiased splicing quantification in the software AltAnalyze (Python2) which accurately identities and quantifies alternative splicing events
from RNA-Seq aligned BAM files. The second component is the separate prediction of MHC-bound peptides (T-antigen) and altered surface proteins (B-antigen) from the observed splicing
events. Hence, the installation of SNAF requires two steps, including the download of necessary dependencies and gene models.

Step 1: AltAnalyze
--------------------

A docker image can be downloaded from DockerHub and run using one line of code::

    # build the image
    docker pull frankligy123/altanalyze
    # run the container, the below command assume a folder named bam is in your current folder on the host machine
    docker run -v $PWD:/usr/src/app/run -t frankligy123/altanalyze bam

The resultant junction count matrix will be in ``./altanalyze_output/ExpressionInput/counts.original.pruned.txt``, all the directory and subdirectory
will be automatically created.

Alternatively, lots of HPC on university or institution use Singularity instead of docker::

    # pull the image
    singularity pull altanalyze.sif docker://frankligy123/altanalyze
    # run the container, the below command assume a folder named bam is in your current folder on the host machine
    singularity run -B $PWD:/usr/src/app/run altanalyze.sif bam


Step 2: SNAF
--------------

SNAF is a python3 package and has been tested on python>=3.7::

    pip install snaf


Step 3: Reference Dataset
---------------------------

The reference datasets include gene sequences, exon intron coordinates and other information, such as a human membrane protein database. Downloading all of
these files will save significant time compared to resorting to REST API while calling::

    curl -o download.tar.gz http://altanalyze.org/SNAF/download.tar.gz
    tar -xzvf download.tar.gz

Step 4: (Optional) Install netMHCpan4.1 and TMHMM2.0
-------------------------------------------------------

.. note:

    Check the Video tutorial for this step: `Install netMHCpan4.1 and TMHMM2.0 for SNAF <https://www.youtube.com/watch?v=KrAzbR5mRIQ>`_.

By default, SNAF uses MHCflurry to predict which peptides will undergo MHC presentation, however, users can optionally install 
netMHCpan4.1 to be used instead. TMHMM2.0 is used for topology prediction in the B-antigen membrane protein workflow if installed. If not installed, results may be less accurate. 
These tools can be downloaded from the authors source websites to get the latest version 
(`netMHCpan4.1 <https://www.cbs.dtu.dk/service.php?NetMHCpan>`_, `TMHMM2.0 <https://services.healthtech.dtu.dk/service.php?TMHMM-2.0>`_). SNAF will ask you
to provide software_path (where do you install these two softwares) when running corresponding steps, that's how these two softwares will be used in SNAF.




