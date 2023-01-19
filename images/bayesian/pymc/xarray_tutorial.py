import numpy as np
import pandas as pd
import xarray as xr
import pymc as pm
import arviz as ar
import pickle

# DataArray
data = xr.DataArray(np.random.randn(2, 3), dims=('x', 'y'), coords={'x': ['x1', 'x2'],'y':['y1','y2','y3']})
print(data.values)
print(type(data.dims))
print(data.coords['x'].to_index())
print(data.attrs)

# Dataset, thinking of data variables that share same definition of dimention, but can have or not have certain dimension with different dtype
ds = xr.Dataset(dict(foo=data, bar=("x", [1, 2]), baz=np.pi))

# InferenceData -- Trace
with open('pickle_mcmc_trace.p','rb') as f:
    trace = pickle.load(f)
'''
Inference data with groups:
	> posterior
	> sample_stats
	> observed_data
'''

posterior = trace.posterior
'''
<xarray.Dataset>
Dimensions:      (chain: 2, draw: 1000, sigma_dim_0: 100, psi_dim_0: 100,
                  mu_dim_0: 100)
Coordinates:
  * chain        (chain) int64 0 1
  * draw         (draw) int64 0 1 2 3 4 5 6 7 ... 993 994 995 996 997 998 999
  * sigma_dim_0  (sigma_dim_0) int64 0 1 2 3 4 5 6 7 ... 92 93 94 95 96 97 98 99
  * psi_dim_0    (psi_dim_0) int64 0 1 2 3 4 5 6 7 8 ... 92 93 94 95 96 97 98 99
  * mu_dim_0     (mu_dim_0) int64 0 1 2 3 4 5 6 7 8 ... 92 93 94 95 96 97 98 99
Data variables:
    sigma        (chain, draw, sigma_dim_0) float64 1.0 1.0 1.0 ... 1.0 1.0 1.0
    psi          (chain, draw, psi_dim_0) float64 0.9214 0.9914 ... 0.9742
    mu           (chain, draw, mu_dim_0) float64 26.18 24.21 ... 24.26 24.97
Attributes:
    created_at:                 2023-01-19T19:02:34.762986
    arviz_version:              0.14.0
    inference_library:          pymc
    inference_library_version:  5.0.1
    sampling_time:              698.9501910209656
    tuning_steps:               1000
'''



