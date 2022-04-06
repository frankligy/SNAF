Development
=============

This section describes the low-level implementation of the program.

At the low level, SNAF follows the Object Oriented Programming (OOP) strategy and construct a series of Python objects.

NeoJunction Class
--------------------

Instantiate and execute the methods associated with each NeoJunction object::

    nj = snaf.NeoJunction(uid='ENSG00000104938:E2.3-E3.1',count=25,check_gtex=False)
    nj.detect_type()
    nj.retrieve_junction_seq()
    nj.in_silico_translation()
    nj.binding_prediction(hlas=hlas,binding_method='netMHCpan')
    nj.immunogenicity_prediction()
    nj.derive_candidates(stage=3,verbosity=1,contain_uid=False)
    nj.visualize(outdir='.',name='check_nj_visualization.pdf')


JunctionCountMatrixQuery Class
--------------------------------

This is basically a wrapper class that we will be dealing with most of the time for T antigen, as the input is a dataframe with all the splicing junction
and samples, this object provide the user interface to process the output from AltAnalyze.


SurfaceAntigen Class
------------------------

For B surface antigen, we can again instantiate and execute the methods::

    uid = 'ENSG00000185499:E7.2-E9.1'
    sa = surface.SurfaceAntigen(uid,False)
    sa.detect_type()
    sa.retrieve_junction_seq()
    sa.recovery_full_length_protein()
    sa.find_orf()
    sa.orf_check(n_stride=2)
    sa.align_uniprot(tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
    sa.visualize(index=12,fragment=None)