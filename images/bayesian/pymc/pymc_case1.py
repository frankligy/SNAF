#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import arviz as az
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pymc3 as pm
import theano

test_scores = pd.read_csv(pm.get_data("test_scores.csv"), index_col=0)
X = test_scores.dropna().astype(float)
y = X.pop("score")
X -= X.mean()
X /= X.std()
N, D = X.shape

D0 = D/2
with pm.Model(coords={'predictors':X.columns.values}) as test_score_model:
    sigma = pm.HalfNormal('sigma',25)
    tau = pm.HalfStudentT('tau',2,D0/(D-D0)*sigma/np.sqrt(N))
    lam = pm.HalfStudentT('lam',2,dims='predictors')
    c2 = pm.InverseGamma('c2',1,0.1)
    z = pm.Normal('z',0,1,dims='predictors')
    beta = pm.Deterministic('beta', z * tau * lam * at.sqrt(c2 / (c2 + tau**2 * lam**2)), dims="predictors")
    beta0 = pm.Normal('beta0',100,25)
    scores = pm.Normal('scores',beta0+at.dot(X.values,beta),sigma,observed=y.values)

pm.model_to_graphviz(test_score_model)
gv.format = 'pdf'
gv.render(filename='model_graph')

with test_score_model:
    prior_samples = pm.sample_prior_predictive(100)
az.plot_dist(
    test_scores["score"].values,
    kind="hist",
    color="C1",
    hist_kwargs=dict(alpha=0.6),
    label="observed",
)
az.plot_dist(
    prior_samples.prior_predictive["scores"],
    kind="hist",
    hist_kwargs=dict(alpha=0.6),
    label="simulated",
)
plt.xticks(rotation=45)

with test_score_model:
    idata = pm.sample(1000, tune=2000, random_seed=42, target_accept=0.99)

az.plot_trace(idata, var_names=["tau", "sigma", "c2"])
az.plot_energy(idata)
az.plot_forest(idata, var_names=["beta"], combined=True, hdi_prob=0.95, r_hat=True)

