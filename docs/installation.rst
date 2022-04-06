Installation
===============

SNAF in generall consists of two parts, the first part is an established algorithm AltAnalyze which accurately identity and quantify alternative splicing events
from RNA-Seq data. The second part is to predict both MHC-bound peptides (T-antigen) and altered surface proteins (B-antigen) from the obtained splicing
events. Hence, the installation also contains two steps, plus downloading necessary data.

Step1: AltAnalyze
--------------------

A docker image can be downloaded from docker Hub and run using one line of code::

    # build the image
    docker pull frankligy123/AltAnalyze
    # run the container
    docker run -v $PWD:/usr/src/app -t frankligy123/AltAnalyze ./path_to_bam_files


Step2: SNAF
--------------

SNAF is a python3 package and has been tested on python>=3.7::

    pip install snaf


Step3: Reference Dataset
---------------------------

The reference datasets include gene sequence, exon intron coordinates and other information like human membrane protein database. Downloading all of
them to the disk will save a lot of time later compared to resorting to REST API while calling::

    curl -o download.tar.gz http://altanalyze.org/SNAF/download.tar.gz
    tar -xzvf download.tar.gz

Step4: (Optional) Install netMHCpan4.1 and TMHMM2.0
-------------------------------------------------------

SNAF doesn't require the presence of these two software, but it is highly recommended to install these two to achieve optimal performance. These two 
tools are super easy to install. Since users have to sign agreement to download the tools, I have to refer you to their official websites 
(`netMHCpan4.1 <https://www.cbs.dtu.dk/service.php?NetMHCpan>`_, `TMHMM2.0 <https://services.healthtech.dtu.dk/service.php?TMHMM-2.0>`_)




