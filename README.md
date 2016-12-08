# SSrGE procedure

This procedure aims to fit sparse linear models using a binary matrix (n_samples x n_SNV) as features matrix and a gene expression matrix (n_genes x n_samples) as response. The procedure a infer sparse linear model (LASSO by default) for each gene (raw in the second matrix) and keep the non-null inferred coefs.

This procedure can be used as dimension reduction/feature selection or feature ranking. It is based on the scikit-learn library and is easy to re-implement. However, the package allows to parallelize the fitting procedures, implements a cross-validation procedure, eeSNVs and gene rankings and can extract SNV and Gene expressions (normalized) matrices from RNA-seq dataset.

SSrGE can be used in stand-alone to reduce any single-cell SNV matrix (raw:single-cell, col: SNV (binary)), using a single-cell gene expression matrix (raw: gene-expression (float), col:single-cell). However, we have developped two additional modules, included in this package, that can be used to download and process RNA-seq data:
* [download_ncbi_data](https://github.com/lanagarmire/SSrGE/README_download_ncbi_rsa.md): download and extract .sra files from NCBI
* [SNV_calling](https://github.com/lanagarmire/SSrGE/README_SNV_calling.md): align reads/infer SNVs and infer gene expression matrices from .fastq files.


## installation (local)

```bash
git clone https://github.com/lanagarmire/SSrGE.git
cd SSrGE
pip install -r requirements.txt --user
```

## Requirements
* Linux/ Unix (not tested) working environment
* [python 2 (>=2.7)](https://www.python.org/download/releases/2.7.2/)
* Python libraries (automatically installed with the pip install command):
  * Numpy
  * Scipy
  * [Scikit-learn](http://scikit-learn.org/) (version = 0.18)
  * tabulate

## usage
* test SSrGE is functional:
```bash
  python test/test_ssrge.py -v
  ```

* Instantiate and fit SSrGE:

SSrGE should be used as a python package, below are usage example.

```python
from garmire_SSrGE.ssrge import SSrGE
from garmire_SSrGE.examples import create_example_matrix_v1 # create examples matrices


help(SSrGE) # See the different functions and specific variables
help(create_example_matrix_v1)

X, Y, W = create_example_matrix_v1()

ssrge = SSrGE()

ssrge.fit(X, Y)

score_models, score_null_models = ssrge.score(X, Y)

X_r = ssrge.transform(X)

print X_r.shape, X.shape

ranked_feature = ssrge.rank_eeSNVs()

ssrge_ES = SSrGE(model='ElasticNet', alpha=01, l1_ratio=0.5) # Fitting using sklearn ElasticNet instead
ssrge_ES.fit(X, Y)

```

* Rank eeSNVs:

```python
ranked_feature = ssrge.rank_eeSNVs()
```

* Performing cross-validation

```python
from garmire_SSrGE.linear_cross_validation import LinearCrossVal

help(LinearCrossVal)

X, Y, W = create_example_matrix_v1()

cross_val = LinearCrossVal(
model='LASSO',
SNV_mat=X,
GE_mat=Y
)

path = cross_val.regularization_path('alpha',  [0.01, 0.1, 0.2])
```

## Use K top-ranked eeSNVs

Instead of relying on the regularization parameter (alpha), to select the number of eeSNVs. One can specify the `nb_ranked_features` argument to obtain a predefined number of eeSNVs, assuming that nb_ranked_features is lower than the number of eeSNVs obtained with the specified alpha.

```python
ssrge_topk = SSrGE(nb_ranked_features=2)
X_r_2 = ssrge_topk.fit_transform(X, Y)

print X_r_2.shape # (100, 2)

```

## Rank genes using eeSNVs and parse SNV ids

In order to rank genes with eeSNVs, the SSrGE instance must be instantiated with SNV ids and gene ids list.

* the gene id order should correspond to the gene matrix
* a SNV id should be a tuple containing the gene id harboring the given SNV and a user defined SNV id (genome position for example).

```python
gene_id_list_example = ['KRAS', 'HLA-A', 'SPARC']
snv_id_list_example = [('KRAS', 10220), ('KRAS', 10520), ('SPARC', 0220)]


## real example
from garmire_SSrGE.examples import create_example_matrix_v2

X, Y, gene_id_list, snv_id_list = create_example_matrix_v2()

ssrge = SSrGE(
      snv_id_list=snv_id_list,
      gene_id_list=gene_id_list,
      nb_ranked_features=2,
      alpha=0.01)

ssrge.fit(X, Y)

print ssrge.ran_genes()

```

## contact and credentials
* Developer: Olivier Poirion (PhD)
* contact: opoirion@hawaii.edu

## related packages
* [SNV calling pipeline](https://github.com/lanagarmire/SNV_calling)
* [NCBI GEO download](https://github.com/lanagarmire/download_ncbi_sra)