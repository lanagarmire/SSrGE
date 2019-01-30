# SSrGE procedure

This procedure aims to fit sparse linear models using a binary matrix (n_samples x n_SNV) as features matrix and a gene expression matrix (n_genes x n_samples) as response. The procedure infers a sparse linear model (LASSO by default) for each gene (raw in the second matrix) and keeps the non-null inferred coefs.

This procedure can be used as  a dimension reduction/feature selection procedure or as a feature ranking. It is based on the Scikit-Learn library and is easy to re-implement. However, the package allows to parallelize the fitting procedures, implements a cross-validation procedure and performs eeSNVs and gene rankings.

SSrGE can be used as a stand-alone procedure to reduce any SNV matrix (raw:single-cell, col: SNV (binary)), using a gene expression matrix (raw: gene-expression (float), col:single-cell). However, we have developped two additional modules, included in this package, that can be used to download and process RNA-seq data:
* [download_ncbi_data](https://github.com/lanagarmire/SSrGE/blob/master/README_download_ncbi_rsa.md): download and extract .sra files from NCBI
* [SNV_calling](https://github.com/lanagarmire/SSrGE/blob/master/README_snv_calling.md): align reads/infer SNVs and infer gene expression matrices from .fastq files.

Alternatively, we compiled the download, alignment, and SNV calling pipelines into a docker container: opoirion/ssrge (see bellow).


## installation (local)

```bash
git clone https://github.com/lanagarmire/SSrGE.git
cd SSrGE
pip2 install -r requirements.txt --user # python 2.7.X must be used
```

## Requirements
* Linux working environment
* [python 2 (>=2.7)](https://www.python.org/download/releases/2.7.2/)
* Python libraries (automatically installed with the pip install command):
  * Numpy
  * Scipy
  * [Scikit-learn](http://scikit-learn.org/) (version = 0.18)
  * tabulate

## usage
* test SSrGE is functional:
```bash
  python2 test/test_ssrge.py -v
  ```

* Instantiate and fit SSrGE:

SSrGE should be used as a python package, below are usage example.
SSrGE takes as input two matrices (A SNV matrix (n_cells x n_SNVs) and a Gene matrix (n_cells x n_Genes)
In the original study, we encoded X with the following procedure:
if a given snv (s) is present into a given cell (c), then X_c,n = 1
However, any type of encoding or continuous values can be used (For example, one can use X_c,n = 1 for a 1/1 genotype and 0.5 for a 0/1 genotype)

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

* Add CNV matrix:

The fit method can take an additional CNV matrix of shape (n_cells x n_genes), and describing the CNV level for each gene.

```python
from garmire_SSrGE.examples import create_example_matrix_v3

X, Y, C, W = create_example_matrix_v3()

help(ssrge.fit) # see the specific documentation of the fit method
ssrge.fit(X, Y, C)
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

Instead of relying on the regularization parameter (alpha), to select the number of eeSNVs, the `nb_ranked_features` argument can be specified to abotained a fixed  number of eeSNVs (assuming that nb_ranked_features is lower than the number of eeSNVs obtained with the specified alpha).

```python
ssrge_topk = SSrGE(nb_ranked_features=2)
X_r_2 = ssrge_topk.fit_transform(X, Y)

print X_r_2.shape # (100, 2)

```

## Ranking genes using eeSNVs and providing SNV ids

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

print ssrge.rank_genes()

```

## Analyzing a subgroup

Extract specific eeSNVs and impacted genes of a given subgroup. a given eeSNV is specific to a subgroup if it is signficantly more present amongst the cells from the given subgroup:

```python

# Defining as a subgroup the first 6 elements from X
subgroup = ssrge.rank_features_for_a_subgroup([0, 1, 2, 3, 4, 5])

print subgroup.ranked_genes
print subgroup.ranked_eeSNVs

print subgroup.significant_genes
print subgroup.significant_eeSNVs

```

## create SNV and GE matrices from .VCF files and gene expression files

It is possible to create an SNV matrix using preexisting .vcf files and also a Gene expression matrix using expression files.

Each cell must have a distinct .vcf file with a unique name (e.g. snv_filtered.vcf) inside a unique folder, specific of the cell, with the name of the cells:

* example:

```bash

data
|-- GSM2259781__SRX1999927__SRR3999457
|   |-- snv_filtered.vcf
|   `-- stdout.log
`-- GSM2259782__SRX1999928__SRR3999458
    |-- snv_filtered.vcf
    `-- stdout.log

```

(stdout.log is not used and were created by the previous analysis)

and similarly for the gene expression files (matrix_counts.txt):

```bash

STAR
|-- GSM2259781__SRX1999927__SRR3999457
|   |-- matrix_counts.txt
|   `-- matrix_counts.txt.summary
`-- GSM2259782__SRX1999928__SRR3999458
    |-- matrix_counts.txt
    `-- matrix_counts.txt.summary

```

(matrix_counts.txt.summary is not used and were created by the previous analysis)

* The format of the expression files supported is the following:

```bash

#gene_name    chromsomes    starting position    ending position    additionnal columns    gene expression
MIR6859-3    chr1;chr15;chr16    17369;102513727;67052   17436;102513794;67119   ...    200
```

* variables (paths and file names) specific to GE and SNV matrix extraction can be defined in the config file: garmire_SSrGE/config.py
* First, a GTF index must be created:

```bash
python2 ./garmire_SSrGE/generate_refgenome_index.py
```

* Once the index generated, the matrices can be genereated easily:


# SRA project download, STAR alignment and SNV calling from scratch using docker

## Requirements
* docker
* possible root access
* 13.8 GB of free memory (docker image) + memory for STAR indexes (usually 20 GB per index) and downloaded data

## installation (local)

```bash
docker pull opoirion/ssrge
mkdir /<Results data folder>/
cd /<Results data folder>/
PATHDATA=`pwd`
```

## usage

The pipeline consists of 3 steps (for downloading the data) and 4 steps for aligning and calling SNVs:

```bash
# Download
docker run --rm opoirion/ssrge download_soft_file -h
docker run --rm opoirion/ssrge download_sra -h
docker run --rm opoirion/ssrge extract_sra -h
# align and SNV calling
docker run --rm opoirion/ssrge star_index -h
docker run --rm opoirion/ssrge process_star -h
docker run --rm opoirion/ssrge feature_counts -h
docker run --rm opoirion/ssrge process_snv -h

```

## example

Let's download and process 2 samples from GSE79457 in a project name test_n2

```bash
# download of the soft file containing the metadata for GSE79457
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge download_soft_file -project_name test_n2 -soft_id GSE79457
# download sra files
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge download_sra -project_name test_n2 -max_nb_samples 2
# exctract sra files
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge extract_sra -project_name test_n2
# rm sra files (optionnal)
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge rm_sra -project_name test_n2
## all these data can also be obtained using other alternative workflows
# here you need to precise which read length to use for creating a STAR index and which ref organism (MOUSE/HUMAN)
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge star_index -project_name test_n2 -read_length 100 -cell_type HUMAN
# STAR alignment
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge process_star -project_name test_n2 -read_length 100 -cell_type HUMAN
# sample-> gene count matrix
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge feature_counts -project_name test_n2
#SNV inference
docker run --rm -v $PATHDATA:/data/results/:Z opoirion/ssrge process_snv -project_name test_n2 -cell_type HUMAN
```


## contact and credentials
* Developer: Olivier Poirion (PhD)
* contact: opoirion@hawaii.edu
