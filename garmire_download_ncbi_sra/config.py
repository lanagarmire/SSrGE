"""
NCBI RNAseq dataset downloader

**** CONFIG FILE ****

"""

from os.path import split as pathsplit
from os.path import abspath

############ Variables #################################################
PROJECT_NAME = "jones_pancreatic_cancer"
# The name of the project (defining the name of the folder)
PATH_DATA = "./example/"
# The absolute path where the project will be created
# and the SRA files downloaded and extracted
PATH_SOFT =  "{0}/{1}.soft".format(PATH_DATA, PROJECT_NAME)
# path toward the .soft file (with the corresponding ftp addresses for the .sra files)
NB_THREADS = 4
# number of threads to use for downloading rsa files
FASTQ_DUMP = "fastq-dump"
# path to the fastq software
FASTQ_DUMP_OPTION = "--split-3 -B"
# options to use to extract the sra (using fastqdump)
# "--split-3 -B is the default" and it is strongly recommended to keep it
LIMIT = None
# define the maximum number of sra files to be downloaded
########################################################################

from garmire_download_ncbi_sra.argv import *
