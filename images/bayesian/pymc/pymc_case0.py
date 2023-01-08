#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import arviz as az
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


## Generate data
alpha = 1
sigma = 1
beta = [1,2.5]
size = 100
rng = np.random.default_rng(8927)
X1 = np.random.randn(size)
X2 = np.random.randn(size) * 0.2
Y = alpha + beta[0] * X1 + beta[1] * X2 + rng.normal(size=size) * sigma

# inference
import pymc3 as pm
basic_model = pm.Model()
with basic_model:
    alpha = pm.Normal('alpha',mu=0,sigma=10)  # stochastic
    beta = pm.Normal('beta',mu=0,sigma=10,shape=2)
    sigma = pm.HalfNormal('sigma',sigma=1)
    mu = alpha + beta[0]*X1 + beta[1]*X2  # determisnitic
    Y_obs = pm.Normal('Y_obs',mu=mu,sigma=sigma,observed=Y)  # observed stochastic
    idata = pm.sample()  # InferenceData
az.plot_trace(idata,combined=True)
plt.savefig('test.pdf',bbox_inches='tight');plt.close()
az.summary(idata,round_to=2).to_csv('test.txt',sep='\t')