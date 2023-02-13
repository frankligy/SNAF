#!/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/bayesian/pytorch_pyro_mamba_env/bin/python3.7

import anndata as ad  # need to install from -c bioconda, not -c conda-forge, but still fail (install 0.6.2 instead), finally use pip to solve (0.8.0)
import numpy as np
import pandas as pd
import os,sys
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
import torch
import pyro
import pickle
import pyro.poutine as poutine

from torch.distributions import constraints
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO


data = torch.zeros(10)
data[0:6] = 1.0

def original_model(data):
    f = pyro.sample("latent_fairness", dist.Beta(10.0, 10.0))
    with pyro.plate("data", data.size(0)):
        pyro.sample("obs", dist.Bernoulli(f), obs=data)

def train(model, guide, lr=0.005, n_steps=201):
    pyro.clear_param_store()
    adam_params = {"lr": lr}
    adam = pyro.optim.Adam(adam_params)
    svi = SVI(model, guide, adam, loss=Trace_ELBO())

    for step in range(n_steps):
        loss = svi.step(data)
        if step % 50 == 0:
            print('[iter {}]  loss: {:.4f}'.format(step, loss))

def model_mle(data):
    # note that we need to include the interval constraint;
    # in original_model() this constraint appears implicitly in
    # the support of the Beta distribution.
    f = pyro.param("latent_fairness", torch.tensor(0.5),
                   constraint=constraints.unit_interval)
    with pyro.plate("data", data.size(0)):
        pyro.sample("obs", dist.Bernoulli(f), obs=data)

def guide_mle(data):
    pass

train(model_mle, guide_mle)
mle_estimate = pyro.param("latent_fairness").item()
print("Our MLE estimate of the latent fairness is {:.3f}".format(mle_estimate))


def guide_map(data):
    f_map = pyro.param("f_map", torch.tensor(0.5),
                       constraint=constraints.unit_interval)
    pyro.sample("latent_fairness", dist.Delta(f_map))

train(original_model, guide_map)
map_estimate = pyro.param("f_map").item()
print("Our MAP estimate of the latent fairness is {:.3f}".format(map_estimate))

