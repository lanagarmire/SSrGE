"""
Instanciate argument passed
Also, helper to print when -H or -h is passed as argument
"""

from sys import argv
from sys import exit

HELP = """
-PROJECT_NAME\tThe name of the project (defining the name of the folder)
-PATH_DATA\tThe absolute path where the project will be created and the SRA files downloaded and extracted
-PATH_SOFT\tpath toward the .soft file (with the corresponding ftp addresses for the .sra files)
-NB_THREADS\tnumber of threads to use for downloading rsa files
-FASTQ_DUMP\tpath to the fastq-dump software
-FASTQ_DUMP_OPTION\toptions to use to extract the sra (using fastq-dump) "--split-3 -B is the default" and it is strongly recommended to keep it
-LIMIT\tdefine the maximum number of sra files to be downloaded (default None)
"""
print('ok')

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

        try:
            exec('{0} = {1}'.format(arg[1:], value))
        except Exception:
            pass
