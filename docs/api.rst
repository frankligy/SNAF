API
=====

API reference manual

MHC bound peptides (T antigen)
---------------------------------

JunctionCountMatrixQuery class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: snaf.snaf.JunctionCountMatrixQuery
    :members: get_membrane_tuples, run, generate_results, visualize

mutation_analysis
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.mutation_analysis

survival_analysis
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.survival_analysis

analyze_neoantigens
~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.analyze_neoantigens

run_dash_T_antigen
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.dash_app.app.run_dash_T_antigen

Surface Antigen (B antigen)
-------------------------------

run
~~~~~~
.. autofunction:: snaf.surface.main.run

generate_results
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.surface.main.generate_results

run_dash_B_antigen
~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.surface.main.run_dash_B_antigen


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

Add gene symbol
~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.downstream.add_gene_symbol_frequency_table


Add chromsome coordinate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: snaf.snaf.add_coord_frequency_table






