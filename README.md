# SSrGE procedure

This procedure aims to fit sparse linear models using a binary matrix (n_samples x n_SNV) as features matrix and a gene expression matrix (n_genes x n_samples) as response. The procedure a infer sparse linear model (LASSO by default) for each gene (raw in the second matrix) and keep the non-null inferred coefs.

This procedure can be used as dimension reduction/feature selection or feature ranking. It is based on the scikit-learn library and is easy to re-implement. However, the package allows to parallelize the fitting procedures, implements a cross-validation procedure, eeSNVs and gene rankings and can extract SNV and Gene expressions (normalized) matrices from RNA-seq dataset.

This package can be used in stand-alone to reduce any single-cell SNV matrix (raw:single-cell, col: SNV (binary)), using a single-cell gene expression matrix (raw: gene-expression (float), col:single-cell). However, we have developped two additional packages that can be used to download and process RNA-seq data:
* [download_ncbi_data](https://github.com/lanagarmire/download_ncbi_sra): download and extract .sra files from NCBI
* [SNV_calling](https://github.com/lanagarmire/SNV_calling): align reads/infer SNVs and infer gene expression matrices from .fastq files.


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
  * [Scikit-learn](http://scikit-learn.org/)
  * tabulate

* To extract SNV and Gene expression matrices from RNA-seq dataset, it is required to:
  * define the corresponding .gtf file of the project (reference gene annotation)
  * For each sample:
    * have raw gene expression inferred by featureCounts, in a single file present in a distinct folder, named according to the sample Id.
    * have VCF file, corresponding to the inferred SNV in a single file present in a distinct folder, named according to the sample Id.
    * The first fields of the .vcf files should correspond to the following example:

      ```text
      chrID     start   SNVid       original    new           score     validSNV
      chrM    12883   rs23245       C       T       12122.77        PASS
      ```

  * all the folders containing the gene expression matrices must be in a distinct folder
  * all the folders containing the VCF files must be in a distinct folder

* The data extraction procedure is particularly well suited for data produced using [garmire_SNV_calling](https://github.com/lanagarmire/SNV_calling) package

## configuration
All the project variables can be defined into the config file (./garmire_SSrGE/config.py). Also, when using directly python class instances, one could access to variables and functions description using the interactive help (see usage) with ipython.


## usage
* test SSrGE is functional:
```bash
  python test/test_ssrge.py -v
  nosetests -v test/test_ssrge.py # alternative using nose
  pytest test/test_ssrge.py -v # alternative using pytest
  ```

* Instantiate and fit SSrGE:

```python
from garmire_SSrGE.ssrge import SSrGE
from garmire_SSrGE.examples import create_example_matrix_v1 # create examples matrices


help(SSrGE) # See the different functions and specific variables
help(create_example_matrix_v1)

X, Y, W = create_example_matrix_v1()

procedure = SSrGE()

procedure.fit(X, Y)

score_models, score_null_models = procedure.score(X, Y)

X_r = procedure.transform(X)

print X_r.shape, X.shape

ranked_feature = procedure.rank_eeSNVs()

procedure_ES = SSrGE(model='ElasticNet', alpha=01, l1_ratio=0.5) # Fitting using sklearn ElasticNet instead
procedure_ES.fit(X, Y)

```

* Rank eeSNVs:

```python
ranked_feature = procedure.rank_eeSNVs()
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

## Extract SNV and GE matrices from RNA-seq dataset:
   Starting from .sra files, one can:
   * first download .sra files directly from NCBI .soft files using our custom package : [download_ncbi_sra](https://github.com/lanagarmire/download_ncbi_sra.git)
   * Then, infer eeSNVs using our SNV calling package : [SNV_calling](https://github.com/lanagarmire/SNV_calling.git)
   * Finally, define the folder locations of the different required inputs:
   ** .soft file
   ** .gtf file
   ** expression matrices folder (same format than output of FeatureCount softwaqre)
   ** vcf folder
   *once all the variables of the project are defined* into the config file (config.py), perform the test:

```bash
  python test/test_dataset.py -v
  nosetests -v test/test_dataset.py # alternative using nose
  pytest test/test_dataset.py -v # alternative using pytest
  ```

* extract matrices

```python
from garmire_SSrGE.extract_matrices_from_dataset import ExtractMatrix

help(ExtractMatrix)

extract_matrix = ExtractMatrix()

SNV_mat = extract_matrix.extract_SNV_mat()
GE_mat = extract_matrix.extract_GE_mat()
```

* perform procedures and rank genes and eeSNVs

```python

procedure.fit(SNV_mat, GE_mat)

ranked_eeSNVS = ssrge.rank_eeSNVs(extract_matrix) # instance of ExtractMatrix is required to obtain eeSNV ids and names

ranked_genes = ssrge.rank_genes(extract_matrix)
```

## contact and credentials
* Developer: Olivier Poirion (PhD)
* contact: opoirion@hawaii.edu