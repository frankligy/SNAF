# SNAF
Splicing Neo Antigen Finder (SNAF) is an easy-to-use Python package to identify splicing-derived tumor neoantigens from RNA sequencing data, it further leverages both deep learning and hierarchical bayesian models to prioritize certain candidates for experimental validations

## Environments

```bash
conda create -n neo_env python=3.7
pip install tensorflow==2.3.0 pandas==1.1.1 numpy==1.18.5
pip install h5py anndata matplotlib seaborn requests xmltodict tqdm
conda install -c conda-forge pymc3 mkl-service   # numpy will be updated to 1.21.3 (seems not, still 1.18.5)
pip install mhcflurry==2.0.5 # for now version seems to not matter
pip install statsmodels lifelines umap plotly   # numba needs to be 0.53
pip install requests xmltramp2 dash-dangerously-set-inner-html
pip install mygene
```
