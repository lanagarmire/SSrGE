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


######## DATA EXTRACTION VARIABLES ############################################
USER = 'opoirion'
#Alias to define the GLOBAL_DATA_ROOT, PROJECT_PATH and PROG_ROOT
# (could be overloaded using reference paths)
PROJECT_NAME = '2015_human_CTC_prostate'
# Project name. Used to create folder
GLOBAL_DATA_ROOT = '/data/{0}/'.format(USER)
# Alias to define the root folder for reference data
# (could be overloaded using reference paths)
PROJECT_PATH = '/home/{0}/data/{1}'.format(USER, PROJECT_NAME)
# Alias to define the output folder
SOFT_PATH = "{0}/{1}/{1}.soft".format(GLOBAL_DATA_ROOT, PROJECT_NAME)
# Absolute path for the .soft file (dataset description) from NCBI
GTF_PATH = '{0}/Illumina_hg19/Annotation/genes.gtf'.format(GLOBAL_DATA_ROOT)
# gtf file of the reference genome
EXPRESSION_MATRIX_FOLDER_PATH = '{0}/expression_profile/STAR/'.format(PROJECT_PATH)
# Path of the folders containing the gene expression matrices
# One distinct folder must be used per sample
GENE_MATRIX_NAME = 'matrix_counts.txt'
# name of the gene expression matrix file
# matrix should be in the same format than the output of FeatureCount software
VCF_FOLDER_PATH = '{0}/snv_pipeline_raw/data/'.format(PROJECT_PATH)
# Path of the folders containing the vcf files
# One distinct folder must be used per sample
VCF_NAME = 'snv_filtered.vcf'
# Name of the VCF the vcf files
##############################################################################
