Change Log
============

Version SNAF.git@e23ce39512a1a7f58c74e59b4b7cedc89248b908 (2024.01.15)
--------------------------------------------------------------------------

I recommend the user who knew our tool from the official publication to use this version, I have tested and it worked. I will not release new pypi version 
in the future, but will keep maintaining the github and release new features using specific commits.


Version 0.7.0 (2023.05.23)
---------------------------

This is a release since the paper's original resubmission, in the past year, we added a ton of new features, main changes are briefly summarized below:

#. differential gene analysis, differential splicing analysis, gene enrichment analysis
#. hla coverage analysis, for each pHLA, what percentage of the population it can cover (benefit)
#. enhanced and more flexible control database, adding two normal cohorts by default, and allowing user to add customized
#. improved algorithms for MLE tumor specificity score, and external BayesTS score (https://github.com/frankligy/BayesTS)
#. improved output for both T and B antigens
#. fix multiple bugs

Version 0.6.0 (2022.07.11)
----------------------------

Adding New functionalities:

#. report B and T antigen in a more readable format
#. allow users to add all info to the frequency table at once
#. add long_read mode for B antigen when inferring full-length isoform and update the B viewer codebase as well


Version 0.5.2 (2022.04.17)
----------------------------

Fixed multiple bugs:

#. scratch folder cleaning issue and multiprocessing conflist issues resolved
#. tumor_specificity issue resolved
#. fix docstring displacement
#. In test folder, contain a complete testing case.

Version 0.5.1 (2022.04.16)
----------------------------

Fixed two bugs:

#. Weckzueg can not be version 2.1.1, force it to be 2.0.2
#. expose add_coord_frequency_table to snaf namespace


Version 0.5.0 (2022.04.06)
----------------------------
Initial release

Prior to Version 0.5.0
--------------------------
I started this project since ``Feb 2020`` during my second rotation in Dr. Nathan Salomonis Lab, The initial code base can be seen from the first GitHub
repository listed below. The neoantigen has been used in a few manuscript and ongoing project, in which I rename it as ``AltNeo-BT``, the codebase can be
seen from the second GitHub repository. 

* `NeoAntigenWorkflow <https://github.com/frankligy/NeoAntigenWorkflow/commit/8aa37114b47513496e0fe14f15155f2bdd159e5d>`_
* `AltNeo-BT <https://github.com/frankligy/AltNeo-BT>`_

