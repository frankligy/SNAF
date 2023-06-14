# Using Kallisto to get gene and transcript TPM

## Build the index

1. Download the transcript fasta from [Gencode v36](https://www.gencodegenes.org/human/release_36.html). This reference is good because
each ENST has its ENSG ID as well, which makes it easier for downstream analysis

2. build the index

```bash
module load kallisto/0.44.0
kallisto index -i transcripts.index transcripts.idx gencode.v36.pc_transcripts.fa
```

3. quantify

```bash
# single end
# bootstrap 100 is the default for EM multi-mapped reads
# l is the fragment size, insert + adaptor, can not know afterwards
# s is the std of fragment size
kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 file.fastq.gz
```

4. understand the results

Here, the effective length is the `l_transcript - l_fragment`, which estimate the position where the fragments can come from a transcript.

See below:

[1] https://www.biostars.org/p/9543067/

[2] https://www.biostars.org/p/252823/


5. summarize to gene level

```python
import os
import sys
import pandas as pd
import numpy as np
import argparse

def main(args):
    path = args.input
    p_dir = os.path.dirname(path)
    df = pd.read_csv(path,sep='\t',index_col=0)
    col1 = []
    col2 = []
    for item in df.index:
        col1.append(item.split('|')[0].split('.')[0])
        col2.append(item.split('|')[1].split('.')[0])
    df['ENST'] = col1
    df['ENSG'] = col2
    df.groupby(by='ENSG')['tpm'].apply(lambda x:x.sum()).to_csv(os.path.join(p_dir,'abundance_gene.tsv'),sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert transcript TPM to gene TPM for kallisto output')
    parser.add_argument('--input',type=str,default=None,help='path to the abundance.tsv file')
    args = parser.parse_args()
    main(args)
```
