# SSrGE procedure

This procedure aims to fit sparse linear models using a binary matrix (n_samples x n_SNV) as features matrix and a gene expression matrix (n_genes x n_samples) as response. The procedure a infer sparse linear model (LASSO by default) for each gene (raw in the second matrix) and keep the non-null inferred coefs.

This procedure can be used as dimension reduction/feature selection or feature ranking. It is based on the scikit-learn library and is easy to re-implement. However, the package allows to parallelize the fitting procedures, implements a cross-validation procedure, vSNVs and gene rankings and can extract SNV and Gene expressions (normalized) matrices from RNA-seq dataset.


## installation (local)

```bash
git clone garmire_SSrGE
cd garmire_SSrGE
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

* The data extraction procedure is particularly well suited for data produced using [garmire_SNV_calling](https://gitlab.com/Grouumf/garmire_SNV_calling) package

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

procedure = SSrGE()

procedure.fit(X, Y)
X_r = procedure.transform(X)

print X_r.shape, X.shape

ranked_feature = procedure.rank_vSNVs()
```

* Rank vSNVs:

```python
ranked_feature = procedure.rank_vSNVs()
```

* Performing cross-validation

```python
from garmire_SSrGE.linear_cross_validation import LinearCrossVal


X, Y, W = create_example_matrix_v1()

cross_val = LinearCrossVal(
model='LASSO',
SNV_mat=X,
GE_mat=Y
)

path = cross_val.regularization_path('alpha',  [0.01, 0.1, 0.2])
```

* Extract SNV and GE matrices from RNA-seq dataset:
  * *once all the variables of the project are defined* into the config file (config.py), perform the test:

```bash
  python test/test_dataset.py -v
  nosetests -v test/test_dataset.py # alternative using nose
  pytest test/test_dataset.py -v # alternative using pytest
  ```