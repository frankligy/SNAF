Tutorial
==========

In this tutorial, we are going to identify the splicing derived neoantigen (both MHC bound T antigen and altered surface B antigen) on 
TCGA Skin Cutaneous Melanoma (SKCM) cohort (472 samples).

.. note::

    I use the whole dataset (472 bam files) to demostrate the full functionalities of SNAF. Actually finishing running this tutorial could take 
    hours and requires multi-core High Performance Computers (HPC). So please replace the bam file folder with your sample and configure your sample 
    HLA type file as well as illustrated below. I hope you find the tutorial intuitive and useful.

Running AltAnalyze to identify alternative splicing event
-----------------------------------------------------------

The journey starts with a bunch of bam files correpoinding to each patient sample, imagine you have all the bam files stored at ``/user/ligk2e/bam``
folder, the full path to the folder is the only thing you need to run this step::

    docker run -v $PWD:/usr/src/app -t frankligy123/AltAnalyze /user/ligk2e/bam

The output of these step contains various useful information for the junction and gene expression, but the most important output that will be used
in the following step is the junction count matrix. Let us take a look at a subsampled junction count matrix, each row represent a splicing junction
annotated by AltAnalyze gene model, and each column represents the sample name. The numerical value represents the number of reads that support the 
occurence of a certain junction.

.. csv-table:: junction count matrix
    :file: ./_static/sample.csv
    :widths: 10,10,10,10,10,10,10,10,10,10,10
    :header-rows: 1

Identify MHC-bound neoantigen (T-antigen)
---------------------------------------------

Now having the junction count matrix, we can go ahead and predict out the MHC-bound neoantigen. The only additional input we need is
the patient HLA type information, in this analyis, we utilized Optitype to infer the 4 digit HLA type from RNA-Seq data, the ``sample_hla.txt`` file 
look like below::

                sample	                                                hla
    TCGA-3N-A9WB-06A-11R-A38C-07.bed	HLA-A*02:01,HLA-A*02:01,HLA-B*39:10,HLA-B*15:01,HLA-C*03:03,HLA-C*12:03
    TCGA-3N-A9WC-06A-11R-A38C-07.bed	HLA-A*02:01,HLA-A*01:01,HLA-B*40:01,HLA-B*52:01,HLA-C*03:04,HLA-C*12:02
    TCGA-3N-A9WD-06A-11R-A38C-07.bed	HLA-A*11:01,HLA-A*32:01,HLA-B*40:02,HLA-B*35:01,HLA-C*04:01,HLA-C*02:02
    TCGA-BF-A1PU-01A-11R-A18S-07.bed	HLA-A*02:01,HLA-A*01:01,HLA-B*07:02,HLA-B*18:01,HLA-C*07:01,HLA-C*07:02


Loading and instantiating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load the packages::

    import os,sys
    import pandas as pd
    import numpy as np
    import snaf

The first step is to load our downloaded reference data into the memory to facilitate the repeat retrieval of the data::

    # database directory (where you extract the reference tarball file)
    db_dir = '/user/ligk2e/download'  
    # instantiate (if using netMHCpan)
    netMHCpan_path = '/user/ligk2e/netMHCpan-4.1/netMHCpan'
    snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
    # instantiate (if not using netMHCpan)
    snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='MHCflurry',software_path=None)

Running T antigen program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first instantiate ``JunctionCountMatrixQuery`` object, here the ``df`` is the junction count matrix (a pandas dataframe) that we showed above.::

    jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df)

We parse the HLA type ``sample_hla.txt`` file a little bit into a nested list, the order of the list should be the same as the samples present in the 
dataframe columns::

    sample_to_hla = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
    hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]

Main program can be wrapped in one line of code, a folder named ``result`` will be created and the resultant ``JunctionCountMatrixQuery``
object will be saved as a pickle file::

    jcmq.run(hlas=hlas,outdir='./result')

To generate a series of useful output including neoantigen burden and neoantigen frequency, we deserialize the pickle file back to memory and automatically
generate these output files::

    snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')

Now in the ``result`` folder, we can have neoantigen burden files associated with each stage, stage correspoindance is as below:

* stage0: neojunction, the amount of tumor-specific junctions
* stage1: peptides that predicted (3-way in-silico translation) from neojunction
* stage2: peptides that are present on MHC molecule (based on netMHCpan or MHCflurry prediction)
* stage3: peptides that are immunogenic (DeepImmuno prediction)

The burden matrix looks like below, the last column and last row represents the mean burden for each feature and total burden for each sample. Since I 
only take the last 10 columns and rows, all of the entries are zero, but I want to show you the summary columns and rows.

.. csv-table:: burden matrix
    :file: ./_static/burden_stage2_sample.csv
    :widths: 10,10,10,10,10,10,10,10,10,10,10
    :header-rows: 1

Neoantigen occurence plot shows the distinctive pattern between shared neoantigens (left part) and unique neoantigens (right part).

.. image:: ./_static/neo_freq.png
    :height: 400px
    :width: 500px
    :align: center
    :target: target

Visualization
~~~~~~~~~~~~~~~~~

A very important question is to ask, what splicing event produce a certain neoepitope, we provide a convenient plotting function to achieve that::

    jcmq.visualize(uid='ENSG00000167291:E38.6-E39.1',sample='TCGA-DA-A1I1-06A-12R-A18U-07.bed',outdir='./result')

.. image:: ./_static/t_visual.png
    :height: 400px
    :width: 500px
    :align: center
    :target: target

Survival Analysis
~~~~~~~~~~~~~~~~~~~~~~~

We download the TCGA SKCM survival data from Xena browser, we provide a convenient function to do survival analyis using various stratification criteria::

    survival = pd.read_csv('TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
    burden = pd.read_csv('result/burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
    burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
    snaf.survival_analysis(burden,survival,n=2,stratification_plot='result/stage2_stratify.pdf',survival_plot='result/stage2_survival.pdf')


.. image:: ./_static/survival.png
    :height: 400px
    :width: 600px
    :align: center
    :target: target

Mutation Association Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We download the TCGA SKCM mutation data from Xena browser, we provide a convenient function to calculate all association and plot them::

    mutation = pd.read_csv('TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
    mutation = mutation.loc[mutation['filter']=='PASS',:]
    burden = pd.read_csv('result/burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
    burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
    snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='result/stage3_mutation.txt')
    snaf.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='result/stage3_mutation_CAMKK2.pdf',genes_to_plot=['CAMKK2'])

.. csv-table:: mutation
    :file: ./_static/stage3_mutation_sample.csv
    :widths: 10,10,10,10
    :header-rows: 1

For a specific mutation ``CAMKK2``, which has been reported that the suppression of this gene can increase the ferroptosis efficacy and 
anti-PD1 immunotherapy (`paper link <https://pubmed.ncbi.nlm.nih.gov/34242660/>`_), we showed that patients with mutated ``CAMKK2`` have higher 
neoantigen burden so that can explain why it lead to better immunotherapy efficacy.

.. image:: ./_static/mutation.png
    :height: 400px
    :width: 600px
    :align: center
    :target: target


Interactive neoantigen Viewer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can launch a dash interactive neoantigen viewer to visualize all the neoantigens based on their physiochemical properties and their motif
composition along with the source splicing junction::

    snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=[9,10],fasta=True)
    snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/shared_vs_unique_neoantigen_all.txt')

.. image:: ./_static/t_viewer.png
    :height: 400px
    :width: 500px
    :align: center
    :target: target


Identify altered surface protein (B-antigen)
-----------------------------------------------

As a separate branch, B-antigen pipeline aims to priotize the altered surface protein from abnormal splicing events.

Instantiating B pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We again load some necessary reference data files to RAM::

    surface.initialize(db_dir=db_dir)

Running the program
~~~~~~~~~~~~~~~~~~~~~~~~~

We first obtain the membrane splicing events::

    membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df)

Then we run the B pipeline::

    # if using TMHMM
    surface.run(membrane_tuples,outdir='result',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
    # if not using TMHMM
    surface.run(membrane_tuples,outdir='result',tmhmm=False,software_path=None)

After this step, a pickle file will again be deposited to the ``result`` folder. However, we do want to generate human-readable results::

    # if having gtf file for long-read data
    surface.generate_results(pickle_path='./result/surface_antigen.p',outdir='result',strigency=5,gtf='./SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf') 
    # if not having 
    surface.generate_results(pickle_path='./result/surface_antigen.p',outdir='result',strigency=3,gtf=None)

Different strigency are explanined below:

* strigency1: novel isoform needs to be absent in uniprot database
* strigency2: novel isoform also needs to be a documented protein-coding gene
* strigency3: novel isoform also needs to not be subjected to Nonsense Mediated Decay (NMD)
* strigency4: novel isoform also needs to have long-read or EST support (as long as the novel junction present in full-length)
* strigency5: novel isoform also needs to have long-read or EST support (whole ORF needs to be the same as full-length)


Interactive neoantigen viewer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to T antigen, users can explore all the altered surface protein for B antigen::

    surface.run_dash_B_antigen(pkl='result/surface_antigen.p',candidates='result/candidates_5.txt',
                               python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')

.. image:: ./_static/b_viewer.png
    :height: 400px
    :width: 600px
    :align: center
    :target: target

Tumor Specificity (GTEx)
----------------------------

For a specific splicing event, we can visualize its tumor specificity by comparing its expression in tumor versus normal tissue::

    snaf.gtex_visual_combine('ENSG00000167291:E38.6-E39.1',norm=True,outdir='result',tumor=df)

here ``norm`` argument controls whether we'd like to normalize the raw read count to Count Per Million (CPM) to account for sequencing depth bias.

.. image:: ./_static/gtex_combine.png
    :height: 400px
    :width: 500px
    :align: center
    :target: target

You can also view each tissue type separately::

    snaf.gtex_visual_subplots('ENSG00000198053:E7.2-E13.1_1915159',norm=True,outdir='result')

.. image:: ./_static/gtex_subplots.png
    :height: 400px
    :width: 500px
    :align: center
    :target: target






