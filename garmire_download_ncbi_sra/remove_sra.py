#! /usr/bin/python

"""
remove sra file
"""

from garmire_download_ncbi_sra.config import PATH_DATA
from os import popen


def main():
    rm_sra()

def rm_sra():
    """extract sra file"""
    path_seq = PATH_DATA + '/fastq/'
    popen('rm {0}*.sra'.format(path_seq)).read()


if __name__ == "__main__":
    main()
