#!/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/bayesian/pytorch_pyro_mamba_env/bin/python3.7

import os,sys
import datetime
from functools import partial
import numpy as np
import torch
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO
from pyro.infer.autoguide import AutoNormal
from pyro.optim import ClippedAdam
import matplotlib as mpl
import matplotlib.pyplot as plt


if torch.cuda.is_available():
    print("Using GPU")
    torch.set_default_tensor_type("torch.cuda.FloatTensor")
else:
    print("Using CPU")

dataset = torch.load('nextstrain.data.pt',map_location=torch.tensor(0.0).device)
counts = dataset['counts']

def get_lineage_id(s):
    """Get lineage id from string name"""
    return np.argmax(np.array([s]) == dataset['lineages'])

def get_location_id(s):
    """Get location id from string name"""
    return np.argmax(np.array([s]) == dataset['locations'])

def get_aggregated_counts_from_locations(locations):
    """Get aggregated counts from a list of locations"""
    return sum([dataset['counts'][:, get_location_id(loc)] for loc in locations])

start = dataset["start_date"]
step = datetime.timedelta(days=dataset["time_step_days"])
date_range = np.array([start + step * t for t in range(len(dataset["counts"]))])

northeast_states = ['USA / Massachusetts',
                    'USA / New York',
                    'USA / Connecticut',
                    'USA / New Hampshire',
                    'USA / Vermont',
                    'USA / New Jersey',
                    'USA / Maine',
                    'USA / Rhode Island',
                    'USA / Pennsylvania']

northeast_counts = get_aggregated_counts_from_locations(northeast_states)  # 27*1316


# The Alpha and Delta variants include many PANGO lineages, which we need to aggregate.
Q_lineages = [lin for lin in dataset['lineages'] if lin[:2] == 'Q.']
AY_lineages = [lin for lin in dataset['lineages'] if lin[:3] == 'AY.']

alpha_lineages = ['B.1.1.7'] + Q_lineages
delta_lineages = ['B.1.617.2'] + AY_lineages

alpha_ids = [get_lineage_id(lin) for lin in alpha_lineages]
delta_ids = [get_lineage_id(lin) for lin in delta_lineages]

alpha_counts = northeast_counts[:, alpha_ids].sum(-1)
delta_counts = northeast_counts[:, delta_ids].sum(-1)
alpha_delta_counts = torch.stack([alpha_counts, delta_counts]).T  # 27*2

def basic_model(counts):
    T, L = counts.shape

    # Define plates over lineage and time
    lineage_plate = pyro.plate("lineages", L, dim=-1)
    time_plate = pyro.plate("time", T, dim=-2)

    # Define a growth rate (i.e. slope) and an init (i.e. intercept) for each lineage
    with lineage_plate:
        rate = pyro.sample("rate", dist.Normal(0, 1))
        init = pyro.sample("init", dist.Normal(0, 1))

    # We measure time in units of the SARS-CoV-2 generation time of 5.5 days
    time = torch.arange(float(T)) * dataset["time_step_days"] / 5.5

    # Assume lineages grow linearly in logit space
    logits = init + rate * time[:, None]

    # We use the softmax function (the multivariate generalization of the
    # sigmoid function) to define normalized probabilities from the logits
    probs = torch.softmax(logits, dim=-1)
    assert probs.shape == (T, L)

    # Observe counts via a multinomial likelihood.
    with time_plate:
        pyro.sample(
            "obs",
            dist.Multinomial(probs=probs.unsqueeze(-2), validate_args=False),
            obs=counts.unsqueeze(-2),
        )

def fit_svi(model, lr=0.1, num_steps=1001, log_every=250):
    pyro.clear_param_store()  # clear parameters from previous runs
    pyro.set_rng_seed(20211214)
    if False:
        num_steps = 2

    # Define a mean field guide (i.e. variational distribution)
    guide = AutoNormal(model, init_scale=0.01)
    optim = ClippedAdam({"lr": lr, "lrd": 0.1 ** (1 / num_steps)})
    svi = SVI(model, guide, optim, Trace_ELBO())

    # Train (i.e. do ELBO optimization) for num_steps iterations
    losses = []
    for step in range(num_steps):
        loss = svi.step()
        losses.append(loss)
        if step % log_every == 0:
            print(f"step {step: >4d} loss = {loss:0.6g}")

    # Plot to assess convergence.
    plt.figure(figsize=(6, 3))
    plt.plot(losses)
    plt.xlabel("SVI step", fontsize=18)
    plt.ylabel("ELBO loss", fontsize=18)
    plt.tight_layout()

    return guide

guide = fit_svi(partial(basic_model, alpha_delta_counts[13:]), num_steps=1501)

for k, v in guide.median().items():
    print(k, v.data.cpu().numpy())


# Model lineage proportions in each region as multivariate logistic growth
def regional_model(counts):
    T, R, L = counts.shape

    # Now we also define a region plate in addition to the time/lineage plates
    lineage_plate = pyro.plate("lineages", L, dim=-1)
    region_plate = pyro.plate("region", R, dim=-2)
    time_plate = pyro.plate("time", T, dim=-3)

    # We use the same growth rate (i.e. slope) for each region
    with lineage_plate:
        rate = pyro.sample("rate", dist.Normal(0, 1))

    # We allow the init to vary from region to region
    init_scale = pyro.sample("init_scale", dist.LogNormal(0, 2))
    with region_plate, lineage_plate:
        init = pyro.sample("init", dist.Normal(0, init_scale))

    # We measure time in units of the SARS-CoV-2 generation time of 5.5 days
    time = torch.arange(float(T)) * dataset["time_step_days"] / 5.5

    # Instead of using the softmax function we directly use the
    # logits parameterization of the Multinomial distribution
    logits = init + rate * time[:, None, None]

    # Observe sequences via a multinomial likelihood.
    with time_plate, region_plate:
        pyro.sample(
            "obs",
            dist.Multinomial(logits=logits.unsqueeze(-2), validate_args=False),
            obs=counts.unsqueeze(-2),
        )












