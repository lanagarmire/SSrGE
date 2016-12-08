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


def main():
    fastq_dump()

def fastq_dump():
    """extract sra file"""
    path_seq = PATH_DATA + 'fastq/'
    for fil in listdir(path_seq):
        if not fnmatch(fil, '*.sra'):
            continue
        print 'go to extraction for file:', fil
        fil = fil.rsplit('.', 1)[0]

        if not isdir("{0}/{1}".format(path_seq, fil)):
            mkdir("{0}/{1}".format(path_seq, fil))
        popen('{3} {2} -v {0}/{1}.sra -O {0}/{1}/'\
              .format(path_seq, fil, FASTQ_DUMP_OPTION, FASTQ_DUMP)).read()


if __name__ == "__main__":
    main()
