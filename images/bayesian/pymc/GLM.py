import anndata as ad
import numpy as np
import pandas as pd
import os,sys
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

import arviz as az
import pymc as pm
import pytensor.tensor as at

import statsmodels.api as sm

# linear regression
# numpy and scipy can solve it as well

spector_data = sm.datasets.spector.load()
spector_data.exog = sm.add_constant(spector_data.exog,prepend=False)
mod = sm.OLS(spector_data.endog,spector_data.exog)
res = mod.fit()
print(res.summary())


# GLM

data = sm.datasets.scotland.load()
data.exog = sm.add_constant(data.exog)
gamma_model = sm.GLM(data.endog,data.exog,family=sm.families.Gamma())
gamma_results = gamma_model.fit()
print(gamma_results.summary())

