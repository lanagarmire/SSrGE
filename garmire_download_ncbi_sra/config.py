"""
NCBI RNAseq dataset downloader

**** CONFIG FILE ****

"""
from sys import argv
from sys import exit


############ Variables #################################################
# The name of the project (defining the name of the folder)
PROJECT_NAME = "chung_notch_2016"
# The absolute path where the project will be created
PATH_DATA = "/data/opoirion/chung_notch_2016/"
# and the SRA files downloaded and extracted
# path toward the .soft file (with the corresponding ftp addresses for the .sra files)
PATH_SOFT =  "{0}/{1}.soft".format(PATH_DATA, PROJECT_NAME)
# number of threads to use for downloading rsa files
NB_THREADS = 4
# number of CPU to be used to extract the .sra files
NB_CPU = 4
# path to the fastq software
FASTQ_DUMP = "fastq-dump"
# options to use to extract the sra (using fastqdump)
# "--split-3 -B is the default" and it is strongly recommended to keep it
FASTQ_DUMP_OPTION = "--split-3 -B"
# define the maximum number of sra files to be downloaded
LIMIT = None
########################################################################


HELP = """
-PROJECT_NAME\tThe name of the project (defining the name of the folder)
-PATH_DATA\tThe absolute path where the project will be created and the SRA files downloaded and extracted
-PATH_SOFT\tpath toward the .soft file (with the corresponding ftp addresses for the .sra files)
-NB_THREADS\tnumber of threads to use for downloading rsa files
-FASTQ_DUMP\tpath to the fastq-dump software
-FASTQ_DUMP_OPTION\toptions to use to extract the sra (using fastq-dump) "--split-3 -B is the default" and it is strongly recommended to keep it
-LIMIT\tdefine the maximum number of sra files to be downloaded (default None)
"""

if argv[1:]:
    if '-h' in argv or '-H' in argv:
        print(HELP)
        exit(0)

    for arg in filter(lambda x:x[0] == '-', argv):

        try:
            value = argv[argv.index(arg) + 1]
        except IndexError:
            continue

        if value.isdigit():
            value = int(value)
        else:
            if value[0] not in ["'", '"'] and value[-1] not in ["'", '"']:
                value = '"{0}"'.format(value)

        try:
            exec('{0} = {1}'.format(arg[1:], value))
        except Exception as e:
            raise(e)
