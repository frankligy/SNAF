import os,sys
from setuptools import setup

# build long description
base_dir = os.path.dirname(os.path.abspath(__file__))
long_description = '\n\n'.join([open(os.path.join(base_dir,'README.md'),'r').read()])

# requires
requires = [
    'tensorflow==2.3.0',
    'pandas==1.3.4',
    'numpy==1.18.5',
    'numba==0.53.0',
    'mhcflurry==2.0.5',
    'h5py==2.10.0',
    'anndata==0.7.6',
    'seaborn==0.11.2',
    'biopython==1.79',
    'requests==2.26.0',
    'xmltodict==0.12.0',
    'xmltramp2==3.1.1',
    'tqdm==4.62.3',
    'scipy==1.4.1',
    'statsmodels==0.13.1',
    'lifelines==0.26.4',
    'umap-learn==0.5.2',
    'plotly==5.4.0',
    'Werkzeug==2.0.2',
    'dash==2.0.0',
    'dash-dangerously-set-inner-html==0.0.2',
    'mygene==3.2.2'
]

setup(
    name = 'SNAF',
    version = '0.7.0',
    description = 'A Python package to predict, prioritize and visualize splicing derived neoantigens, including MHC-bound peptides (T cell antigen) and altered surface protein (B cell antigen)',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    author = 'Guangyuan(Frank) Li',
    author_email='li2g2@mail.uc.edu',
    maintainer='Guangyuan(Frank) Li',
    maintainer_email = 'li2g2@mail.uc.edu',
    url='https://github.com/frankligy/SNAF',
    project_urls={
        'Documentation':'https://snaf.readthedocs.io',
    },
    packages=['snaf','snaf/deepimmuno','snaf/dash_app','snaf/surface'],
    package_data = {'snaf/deepimmuno':['data/*','models/cnn_model_331_3_7/.data-00000-of-00001','models/cnn_model_331_3_7/.index','models/cnn_model_331_3_7/checkpoint'],
                    'snaf':['mqpar.xml','HLA_Allele_frequency_21_populations.csv']},
    install_requires=requires,
    python_requires='>=3.7',
    classifers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 3',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
    ]
)