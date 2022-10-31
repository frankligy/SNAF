API
=====

API reference manual


MHC bound peptides (T antigen)
---------------------------------

initialize
~~~~~~~~~~~
.. autofunction:: snaf.__init__.initialize

JunctionCountMatrixQuery class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: snaf.snaf.JunctionCountMatrixQuery
    :members: get_membrane_tuples, run, generate_results, visualize, show_neoantigen_as_fasta

mutation_analysis
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.mutation_analysis

survival_analysis
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.survival_analysis

survival_regression_cox_model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.survival_regression

analyze_neoantigens
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.analyze_neoantigens

run_dash_T_antigen
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.dash_app.app.run_dash_T_antigen

.. _reference_to_report_candidates:

report_candidates
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.report_candidates

.. _reference_to_proteomics:

Interface to proteomics
----------------------------

chop_normal_pep_db
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.proteomics.chop_normal_pep_db

compare_two_fasta
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.proteomics.compare_two_fasta

remove_redundant
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.proteomics.remove_redundant

set_maxquant_configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.proteomics.set_maxquant_configuration

Surface Antigen (B antigen)
-------------------------------

initialize
~~~~~~~~~~~~~
.. autofunction:: snaf.surface.main.initialize

run
~~~~~~
.. autofunction:: snaf.surface.main.run

generate_results
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.surface.main.generate_results

run_dash_B_antigen
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.surface.main.run_dash_B_antigen

.. _reference_to_report_B_candidates:

report_candidates
~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.surface.main.report_candidates


GTEx Viewer (tumor specificity)
------------------------------------

gtex_visual_combine
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.gtex_viewer.gtex_visual_combine  

gtex_visual_subplots
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.gtex_viewer.gtex_visual_subplots

Miscellaneous
------------------

.. _reference_to_add_gene_symbol:

Add gene symbol
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.add_gene_symbol_frequency_table

.. _reference_to_add_chromsome_coordinate:

Add chromsome coordinate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.snaf.add_coord_frequency_table

.. _reference_to_add_specificity:

Add tumor specificity scores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.gtex.add_tumor_specificity_frequency_table

Add All three info at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.snaf.enhance_frequency_table

Prepare DEG analysis
~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.prepare_DEG_analysis

Visuliaze DEG result (Volcano)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.visualize_DEG_result

Prepare GO analysis
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.prepare_GO_analysis

Visualize GO result
~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.visualize_GO_result







