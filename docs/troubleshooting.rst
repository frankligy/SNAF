Troubleshooting
================

I will put Frequently asked question here and provide solutions. Please refer to the Github Issue page for all the issues.

1. Singularity Error
----------------------

If you encounter problem like below::

    Current folder is /run
    BAM file folder is bam
    /usr/src/app/AltAnalyze.sh: line 25: cd: bam: No such file or directory

That is because there's a folder called ``/run`` in your singularity sandbox, so the first line of ``AltAnalyze.sh`` which use 
relative path to enter /run folder will mistakenly go into the root folder, in turn cause the bam folder not found. The solution is
go to your altanalyze/ sandbox, in /usr/src/app folder edit the AltAnalyze.sh first line from ``cd /run`` to ``cd /usr/src/app/run``. Now 
the problem should be solved.