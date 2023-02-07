#!/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/bayesian/pytorch_pyro_mamba_env/bin/python3.7

import os,sys
import torch
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyro
import logging

pyro.enable_validation(True)
pyro.set_rng_seed(1)
logging.basicConfig(format='%(message)s', level=logging.INFO)

# I want to generate a plate notation for the bayesian NN
import pyro.distributions as dist
import pyro.distributions.constraints as constraints

sim_rna = torch.randint(low=0,high=100,size=(1000,20000))
batch = torch.tensor(np.random.choice(10,size=1000))
def aim3_model(sim_rna,batch):
    n,g = sim_rna.shape
    t = len(torch.unique(batch))

    # underlying u
    u_ng = pyro.param('u_ng',lambda:torch.randn((n,g)),constraint=constraints.positive)

    # dispersion
    sigma_g_beta = pyro.sample('beta_o',dist.Gamma(concentration=9,rate=3))
    sigma_g = pyro.sample('o_g',dist.Exponential(rate=sigma_g_beta).expand([g]))
    
    # additive shift
    mu_t = pyro.sample('mu_t',dist.Gamma(1,100).expand([t]))
    beta_g = pyro.sample('beta_g',dist.Gamma(concentration=9,rate=3).expand([g]))
    o_tg = pyro.sample('o_tg',dist.Exponential(rate=beta_g).expand([t,g]))
    alpha_tg = 1/o_tg**2
    rate = alpha_tg/mu_t.unsqueeze_(1)
    s_tg = pyro.sample('s_tg',dist.Gamma(concentration=alpha_tg,rate=rate))
    s_ug = s_tg[batch,:]

    w_ng = pyro.sample('w_ng',dist.Gamma(concentration=u_ng+s_ug,rate=1/sigma_g**2))  # automatic broadcase as the rightmost are aligned
    with pyro.plate('obs'):
        rna_samples = pyro.sample('rna',dist.Poisson(rate=w_ng),obs=sim_rna)
    return rna_samples
pyro.render_model(aim3_model, model_args=(sim_rna,batch), render_distributions=True, render_params=True, filename='aim3.pdf')
