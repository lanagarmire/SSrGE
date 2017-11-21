#! /usr/bin/python

"""
remove sra file
"""
from os.path import isdir
from os import mkdir
from os import listdir
from os import popen

from fnmatch import fnmatch

from garmire_download_ncbi_sra.config import PATH_DATA
from garmire_download_ncbi_sra.config import FASTQ_DUMP
from garmire_download_ncbi_sra.config import FASTQ_DUMP_OPTION
from garmire_download_ncbi_sra.config import LIMIT
from garmire_download_ncbi_sra.config import NB_CPU


from multiprocessing import Pool


############ VARIABLE ############
PATH_SEQ = PATH_DATA + 'fastq/'
##################################


def main():
    fastq_dump()

def fastq_dump():
    """extract sra file"""
    count = 0

    print('extracting .sra files into: {0}'.format(PATH_SEQ))

    file_list = []

    for fil in listdir(PATH_SEQ):
        if not fnmatch(fil, '*.sra'):
            continue

        file_list.append(fil)

        count += 1

        if LIMIT and count > LIMIT:
            break

    pool = Pool(NB_CPU)

    pool.map(_fastq_dump, file_list)

def _fastq_dump(fil):
    """ """
    print('go to extraction for file:', fil)
    fil = fil.rsplit('.', 1)[0]

    if not isdir("{0}/{1}".format(PATH_SEQ, fil)):
        mkdir("{0}/{1}".format(PATH_SEQ, fil))
    popen('{3} {2} -v {0}/{1}.sra -O {0}/{1}/'\
          .format(PATH_SEQ, fil, FASTQ_DUMP_OPTION, FASTQ_DUMP)).read()


if __name__ == "__main__":
    main()
