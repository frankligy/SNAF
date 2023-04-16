## splice-inferelator project

1. Building the prior based on comprehensive perbutation data on splicing factors, including two cell lines (K562 and HepG2) across two technologies (CRISPR and RNAi).

We started with downloading the fastq files. The raw metadata are in this folder, we processed them a bit using

```bash
./scripts/test_controls_complete.py
```

to make sure all control sampels are in control tsv and only useful control fastq will be downloaded. the `processed_metadata` will be used in the following analysis.

Now generated the `fastq.txt` file to dowload based on

```bash
./scripts/download.py
```

for downloading (2730 fastq files), now start the actually downloading

```bash
# in sub_all.sh files, make sure to employ 10 cores
cd ../fastq
cat ../fastq.txt | xargs -L 1 -P 10 curl -O -J -L 

```
After downloading, check md5sum, this can be a multiple round process, have to check until right

```bash
./scripts/inspect_md5sum.py
```

Then, using `pair_fastq.py` to build the pair, so we have a `pair_info_fastq.txt` file, 1365 pairs.

Now, running `star_grch38.sh` for alignment, construct a file called `STAR_output`.

After that, using `bam_to_bed.sh` to obtain all the junction bed, move the bed file to a folder called `junction`.

From here, we can start running multipath-psi using `multipath_psi.sh`, by convention, I create a folder called `altanalyze_output`. Here, we use `build_group_and_comp.py` to build group and comp files according to different requirements. For examples, first round we just need a dummy group file because we are not going to really look into the metaDataAnalysis results. We also need to build group file for iterative metaDataAnalysis for both "within" and "combine" mode. Additionally, we build group file for all exp and control group, which will be useful in the batch effect inspection.

We use `any_batch_effect.py` to explore batch effect, using EventAnnotation file. We used limma to remove batch effect, here involves the `batch_correction` folder where the `batch.py` script is to build batch information on each category. We have a `companion.py` python script to accompany limma R code because some operations are slow in R. `run_from_count.sh` will rescue the corrected EventAnnotation file. Looking at the `batch_correct` folder within `altanalyze_output` for the batch effect inspection results.

We run metaDataAnalysis iteratively, by doing so, we comduct all the actual computation in a seperate folder `test`, where we copy and pasted all the correponding EventAnnotation file. We use `metaDataAnalysis.sh` to automate the whole process, the resultant PSI matrix will be in `das_within` and `das_combine` folder in `altanalyze_output`. 

The concordance analysis is conducted using `concordance.py` script. 

2. Now, we first focus on shRNA_K562, build functional prior and binding prior from the KD and eclip

Now we utilized `functional_prior.py` to build the prior network from functional evidence. The prior is in `altanalyze_output/prior/functional_prior`. Using `visualize.py` to visualize the prior network as an interactive graph.

For binding prior, the First step is to download the bed file from ENCODE, we are going to relex the IDR peak a bit, so we download the bed files for each replicate. We have the `bed_url.txt` file, and we added teh double quotations to get `bed_url_double_quoted.txt` for downloading the bed files. 

Building `idr_new_env` as a seperate conda env, download the `idr-2.0.4` package, remember you have to download the 2.0.4, and I can not recall where I find the correct package, but it's on the clutser now. We need to call run `idr` programmatically by pairing the two matched bed files, so we wrote a script `assist.py` to run it, and all the idr_bed will be in `idr_bed` folder for further usage.

The constructed binding prior will be in `../altanalyze_output/prior/binding_prior` folder.

The whole prior was produced by mergeing the two sources of priors, we first try rule-based approach, see `whole_prior.py` script. The constructed whole prior will be in `../altanalyze_output/prior/whole_prior` folder.

3. To validate, we want to apply it on HepG2 shRNA, using whole prior built on top of K562. Because it has largest overlap interms of common SF that underwent KD experiments, and decent amount of overlap for e-clip validation. The validation contains several parts:

- derive TFA
- binding evidence
- functional validation





