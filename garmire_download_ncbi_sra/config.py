"""
NCBI RNAseq dataset downloader

**** CONFIG FILE ****

"""
from argparse import ArgumentParser

ARGPARSER = ArgumentParser(description='Argument for the SRA downloading pipeline',
                                   prefix_chars='-')

ARGPARSER.add_argument('-project_name',
                       help='name of the project folder and where to find the fastq files (default: sample_test)',
                       default="sample_test",
                       metavar='str')

ARGPARSER.add_argument('-dl_nb_threads',
                       help=' number of CPU to be used to extract the .sra files (default: 4)',
                       default=4,
                       type=int,
                       metavar='int')

ARGPARSER.add_argument('-nb_cpus',
                       help=' number of CPU to be used to extract the .sra files',
                       default=4,
                       type=int,
                       metavar='int')

ARGPARSER.add_argument('-max_nb_samples',
                       help=' max number of samples downloaded (default None)',
                       default=0,
                       type=int,
                       metavar='int')

ARGPARSER.add_argument('-soft_id',
                       help=' SRA ID used to download the corresponding .soft file (example: "GSE79457")',
                       default="",
                       type=str,
                       metavar='str')

ARGS = ARGPARSER.parse_known_args()[0]

############ Variables #################################################
# The name of the project (defining the name of the folder)
PROJECT_NAME = ARGS.project_name
# The name of the project (defining the name of the folder)
PROJECT_NAME = ARGS.project_name
# The absolute path where the project will be created
PATH_DATA = "/data/results/{0}".format(PROJECT_NAME)
# and the SRA files downloaded and extracted
# path toward the .soft file (with the corresponding ftp addresses for the .sra files)
PATH_SOFT =  "{0}/{1}.soft".format(PATH_DATA, PROJECT_NAME)
# number of threads to use for downloading rsa files
NB_THREADS = ARGS.dl_nb_threads
# number of CPU to be used to extract the .sra files
NB_CPU = ARGS.nb_cpus
# path to the fastq software
FASTQ_DUMP = "fastq-dump"
# options to use to extract the sra (using fastqdump)
# "--split-3 -B is the default" and it is strongly recommended to keep it
FASTQ_DUMP_OPTION = "--split-3 -B"
# define the maximum number of sra files to be downloaded
LIMIT = ARGS.max_nb_samples
#soft ID
SOFT_ID = ARGS.soft_id
########################################################################
