Installation
===============

SNAF in generall consists of two parts, the first part is an established algorithm AltAnalyze which accurately identity and quantify alternative splicing events
from RNA-Seq data. The second part is to predict both MHC-bound peptides (T-antigen) and altered surface proteins (B-antigen) from the obtained splicing
events. Hence, the installation also contains two steps.

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

