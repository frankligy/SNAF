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


2. AltAnalyze warnings and errors
--------------------------------------

You can ignore messages shown below, it won't affect your program.

Something like that::

    library lxml not supported. WikiPathways and LineageProfiler visualization will not work. Please install with pip install lxml.
    Traceback (most recent call last):
    File "/usr/src/app/altanalyze/import_scripts/BAMtoJunctionBED.py", line 217, in parseJunctionEntries
        try: splicesite_db,chromosomes_found, gene_coord_db = retreiveAllKnownSpliceSites()
    File "/usr/src/app/altanalyze/import_scripts/BAMtoJunctionBED.py", line 133, in retreiveAllKnownSpliceSites
        for file in os.listdir(parent_dir):
    OSError: [Errno 2] No such file or directory: 'S'

Or::

    AltAnalyze depedency not met for: patsy

    ...Combat batch effects correction requires pandas and patsy

    AltAnalyze depedency not met for: wx

    ...The AltAnalyze Results Viewer requires wx

    AltAnalyze depedency not met for: fastcluster

    ...Faster hierarchical cluster not supported without fastcluster

    AltAnalyze depedency not met for: lxml

    ...Wikipathways visualization not supported without lxml

    

    WARNING!!!! Some dependencies are not currently met.

    This may impact AltAnalyze's performance


