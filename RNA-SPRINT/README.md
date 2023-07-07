# RNA-SPRINT

RNA-SPRINT is an approach to estimate the RNA binding protein (RBP) activity in each sample. The idea of RNA-SPRINT is simple, we first built a prior network representing the confident relationship between RBP and its targerted splicing junction in K562 data from large-scale knock-down evidence and eCLIP-Seq binding evidence. Using this prior network as basis, we ran Multivariate Decision Tree (MDT) from `decoupler-py` package, which shows consistent better performance in the independent benchmark in HepG2 cell line inference.

The input of RNA-SPRINT is the EventAnnotation file from AltAnalyze results (n_splicing_junctions * n_samples), the output is a `mdt_estimate` file (n_RBPs * n_samples) containing the activity for each RBP in each sample.

We provided a `yml` file to reproduce the environment, but you can build youself based on the `RNA-SPRINT.py` code, should be pretty straightforward.

To run that, you just need to change the capital variables indicated in the `RNA-SPRINT.py` file.