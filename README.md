# MutModels

description of the mutmodels package

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Installation](#installation)
- [Usage](#usage)
  - [Two-sample DM bootstrap test on trinucleotide frequencies](#two-sample-dm-bootstrap-test-on-trinucleotide-frequencies)
    - [mutmodels.dm\_twosample\_SBS](#mutmodelsdm_twosample_sbs)
  - [One-sample DM bootstrap test on trinucleotide frequencies](#one-sample-dm-bootstrap-test-on-trinucleotide-frequencies)
    - [mutmodels.dm\_onesample\_SBS](#mutmodelsdm_onesample_sbs)
  - [custom (mutation) count data](#custom-mutation-count-data)
- [Citation](#citation)
- [Copyright](#copyright)
- [Contact Information](#contact-information)

## Installation
Clone or download this repository. Then, from the base folder of this repository, run:

```bash
pip install .
```
## Usage

### Two-sample DM bootstrap test on trinucleotide frequencies

#### mutmodels.dm_twosample_SBS
Compare the trinucleotide mutation spectrum of two sets of samples
```python
from mutmodels import dm_twosample_SBS

matrix_fn = "./test_data/fig2.SBS96.all" # table with counts per sample
bg_fn = './test_data/fig2.UNT_background.SBS96.tsv' # background signature to subtract, optional

# groups of samples to compare
g2 = ['JG0139b','JG0139c','JG0139d','JG0139e']
g1 = ['JG0140b','JG0140c','JG0140d','JG0140e']

# select subset of substitutions to include in test, optional
substitutions = ['T>G','T>C','T>A']

result = dm_twosample_SBS(matrix_fn,g1,g2,bg_fn=bg_fn, matrix_type='SBS96',
                          substitutions=substitutions,rng=1234)

print(result['p_value'])
print(result['TVD'])
```
| Parameter | Type | Description |
| --- | --- | --- |
| `matrix_fn` | string | path to .tsv-format matrix of MutationType counts per sample  |
| `g1`,`g2` | list | sample ids per group, corresponding to columns in 
| `bg_fn` | string | path to matrix with 1 column containing expected background signal, to be subtracted from sample counts. default=`None` |
|`substitutions`|list (e.g.`['C>T','T>N']`)|optional subset of substitutions to be used for set. default=None (all substitutions considered).|
|`matrix_type`|string | mutational context type. currently only "SBS96" is supported. default=`None`
|`n_bootstraps`| int | number of parametric bootstrap iterations performed. default=`1000`
|`sig_level`| float | alpha-style significance threshold. default=`0.05`
|`dispersion_type`|`'split'`,`'shared'` or `'highest'`| determines method of assigning concentration parameters. Separate concentration parameter fit for each sample group (`'split'`), single parameter fit on the pooled set of samples used for both (`'shared'`), or separate parameter fir for each group, highest dispersion is used for both groups (`'highest'`).
|`studentize`| `True` or `False`| Whether the metric is absolute difference or moderated studentized difference. default = `True`
|`stat_type`| `'max'`(default) or `'sum'` | max-type statistic of summed statistic. default = `'max'`
|`transform`| function or `None` | transformation to perform on proportions before calculating statistic. default=`None`
|`rng`|`None`,int or np.rng instance| for reproducibility of bootstrapping procedure. default=`None`. 

### One-sample DM bootstrap test on trinucleotide frequencies

#### mutmodels.dm_onesample_SBS
Compare the trinucleotide mutation spectrum of a set of samples to a reference signature (e.g. a cosmic mutational signature)

```python
from mutmodels import dm_onesample_SBS

# table with MutationType counts per sample
matrix_fn = './test_data/fig3.corrected.SBS96.all'

# table with MutationType frequencies of signature
ref_sig_fn = './test_data/cosmic_v3.5_GRCh38_SBS30.tsv'


# no set of samples provided, so all will be used
# only trinucleotide spectrum of C>T mutations will be considered
# no studentized test (default), absolute difference in observed and signature proportions used
result = dm_onesample_SBS(matrix_fn,ref_sig_fn,substitutions=['C>T'],
                          matrix_type='SBS96',studentize=False,rng=1234)

print(result['p_value'])
print(result['TVD'])
```
| Parameter | Type | Description |
| --- | --- | --- |
| `matrix_fn` | string | path to .tsv-format matrix of MutationType counts per sample  |
| `ref_fn` | string | path to matrix with 1 column containing reference frequencies, to compare observed data to. default=`None` |
| `sample_ids` | list | sample ids per group, corresponding to columns in matrix_fn |
| `bg_fn` | string | path to matrix with 1 column containing expected background signal, to be subtracted from sample counts. default=`None` |
|`substitutions`|list (e.g.`['C>T','T>N']`)|optional subset of substitutions to be used for set. default=None (all substitutions considered).|
|`matrix_type`|string | mutational context type. currently only "SBS96" is supported. default=`None`
|`n_bootstraps`| int | number of parametric bootstrap iterations performed. default=`1000`
|`sig_level`| float | alpha-style significance threshold. default=`0.05`
|`studentize`| `True` or `False`| Whether the metric is absolute difference or moderated studentized difference. default = `True`
|`stat_type`| `'max'`(default) or `'sum'` | max-type statistic of summed statistic. default = `'max'`
|`transform`| function or `None` | transformation to perform on proportions before calculating statistic. default=`None`
|`rng`|`None`,int or np.rng instance| for reproducibility of bootstrapping procedure. default=`None`. 

### custom (mutation) count data
In case you have your own mutation classification or sets of counts of any other nature, the two tests can be directly run on vector(s) of counts like so:

```python
import numpy as np
from mutmodels.dirichlet_multinomial.core import dm_onesample, dm_twosample

# insert your count vectors as numpy arrays here
g1_counts = [...]  # rows are samples
g2_counts = [...]  # columns are (mutation) categories 

twosample_result = dm_twosample(
    g1_counts,g2_counts,n_bootstraps=100000,rng=1234
)

print(twosample_result['p_value'])
print(twosample_result['TVD'])

```

```python
import numpy as np
from mutmodels.dirichlet_multinomial.core import dm_onesample, dm_twosample

obs_counts = [...]
ref_proportions = [...]

# summed statistic over categories instead of default max-type statistic
onesample_result = dm_onesample(
    obs_counts,ref_proportions,stat_type='sum'
)

print(onesample_result['p_value'])
print(onesample_result['TVD'])
```

## Citation
Manuscript pending publication

## Copyright
[pending]

## Contact Information

For Any questions, requests or bug reports, please contact Joeri van Strien at: j.strien@hubrecht.eu
