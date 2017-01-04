from garmire_SNV_calling.config import ANNOTATION_PATH
from garmire_SNV_calling.config import PATH_OUTPUT as PROJECT_PATH
from garmire_SNV_calling.config import SOFT_PATH as SNV_CALLING_SOFT_PATH
from garmire_SNV_calling.config import MATRIX_OUTPUT_PATH
from garmire_SNV_calling.config import OUTPUT_PATH_SNV

"""
CONFIG file for SSrGE

Principal default values for SSrGE class parameters

The config file gives also parameters in order to
extract SNV and GE matrices from a given project

"""

######## SSrGE VARIABLE ##############################
TIME_LIMIT = 5
# time limit for one linear regression model
MIN_OBS_FOR_REGRESS = 10
# Min number of cell having non null gene expression
# to infer a sparse linear model
NB_THREADS = 4
# Number of threads to run in parallel
CROSS_VAL_NFOLD = 5
# Number of folds to perform the cross validation
######################################################


##################################### vcf and ecpression data #######################
# Paths used to create the SNVs and the gene expression matrices
# The folder architecture used by default is the one from the SNV calling package
# see the file ./garmire_SSrGE/garmire_SNV_calling/config.py
# All the paths defined bellow can be overwritted using a user defined path instead

EXPRESSION_MATRIX_FOLDER_PATH = '{0}/STAR/'.format(MATRIX_OUTPUT_PATH)
# the path for the folders containing the expression matrix files
# one folder per single cell and each folder contains a unique expression matrix (.txt) file
GENE_MATRIX_NAME = 'matrix_counts.txt'
# the name of the file containing the expression matrix inside single-cell each folder
# one folder per single cell and each single-cell folder contains a unique expression file
VCF_FOLDER_PATH = '{0}/data/'.format(OUTPUT_PATH_SNV)
# the path for the folders containing the .vcf files
# one folder per single cell and each single-cell folder contains a unique .vcf file
VCF_NAME = 'snv_filtered.vcf'
# the name of the file containing the vcf inside each folder
INDEX_SAVE_PATH = "{0}/gtf_index/".format(PROJECT_PATH)
# path to save the GTF index
GTF_PATH = ANNOTATION_PATH
# path used in the SNV_calling module
SOFT_PATH = SNV_CALLING_SOFT_PATH # OPTIONNAL, path of the .soft file from ncbi
######################################################################################
