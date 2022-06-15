[![Documentation Status](https://readthedocs.org/projects/snaf/badge/?version=latest)](https://snaf.readthedocs.io/en/latest/?badge=latest)  [![Pypi](https://img.shields.io/pypi/v/snaf?logo=PyPI)](https://pypi.org/project/snaf/)  [![Downloads](https://pepy.tech/badge/snaf)](https://pypi.org/project/snaf/)  [![Stars](https://img.shields.io/github/stars/frankligy/SNAF)](https://github.com/frankligy/SNAF/stargazers)

# SNAF
Splicing Neo Antigen Finder (SNAF) is an easy-to-use Python package to identify splicing-derived tumor neoantigens from RNA sequencing data, it can 
predict, prioritize and visualize MHC-bound neoantigen for T cell (T antigen) and altered surface protein for B cell (B antigen).

![workflow](./images/fig1.png)


# Tutorial and documentation

[Full Documentation](https://snaf.readthedocs.io)

# Input and Output

Simply put, user needs to supply ``a folder with bam files``, and the ``HLA type`` assciated with each patient (using your favorite HLA typing tool). And it will generate predicted immunogenic MHC-bound peptides and altered surface protein. Moreover, there's a myriad of convenient function that enables users to conduct survival analysis, association analysis and publication-quality visualiztion. Check our tutorials for more detail.

# Special notes 

This is currently alpha version pre-release for AACR conference, I am trying to glean feedbacks from AACR attendees and improve my pipeline. If you have any confusion about my desription in the tutorial, run into troubleS while using it, Please feel free to reach out to me and I will be responsive. We are also open to any form of colloborations!

AACR Poster can be found [here](./images/AACR-Li_56x36_CCH_UC_White.pdf).

# Citation

[Guangyuan Li, Nathan Salomonis. SNAF: Accurate and compatible computational framework for identifying splicing derived neoantigens [abstract]. In: Proceedings of the American Association for Cancer Research Annual Meeting 2022; 2022 Apr 8-13. Philadelphia (PA): AACR; Cancer Res 2022;82(12_Suppl):Abstract nr 1898.](https://aacrjournals.org/cancerres/article/82/12_Supplement/1898/701846/Abstract-1898-SNAF-Accurate-and-compatible)

A preprint will be released soon.

# Contact

Guangyuan(Frank) Li

Email: li2g2@mail.uc.edu

PhD student, Biomedical Informatics

Cincinnati Children’s Hospital Medical Center(CCHMC)

University of Cincinnati, College of Medicine