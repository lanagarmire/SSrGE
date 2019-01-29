#! /usr/bin/python

"""
download data from NCBI according to GEO accession number
"""

from os.path import isdir
from os.path import isfile
from os import popen

import urllib2
import re

from distutils.dir_util import mkpath

from garmire_download_ncbi_sra.config import PATH_SOFT
from garmire_download_ncbi_sra.config import PATH_DATA
from garmire_download_ncbi_sra.config import NB_THREADS
from garmire_download_ncbi_sra.config import LIMIT

from multiprocessing.pool import ThreadPool

from time import sleep


############ VARIABLES ############
PATH_SEQ = PATH_DATA + '/fastq/'
###################################


def main():
    download_data()

def _download(url):
    """ """
    gsm, address = url

    try:
        srx = address.rsplit('/', 1)[-1]
        url = urllib2.urlopen(address).read().split()
        srr = url[-1]

        srr_url = "{0}/{1}/{1}.sra".format(address, srr)
        f_name = "{0}{1}__{2}__{3}.sra".format(PATH_SEQ,
                                               gsm,
                                               srx,
                                               srr)
    except Exception as e:
        print('error with SRX {0}!!!'.format(address))
        return "{1} {0}".format(str(e), address)

    if isfile(f_name):
        print("{0} already exists continue...".format(f_name))
        return "{0} already exists continue...".format(f_name)

    try:
        print('downloading {0} to {1}...'.format(srr_url, f_name))
        popen("wget -O {0} {1} --no-verbose".format(
            f_name,
            srr_url)).read()
        print('{0} successfully downloaded'.format(f_name))
        return

    except Exception as e:
        print('error while downloading {0}!!!'.format(address))
        return "{1} {0}\n".format(str(e), address)

    sleep(0.2)

def download_data():
    """download dataset from ncbi """

    urls = get_urls()

    if LIMIT:
        urls = urls[:LIMIT]

    if not isdir(PATH_SEQ):
        mkpath(PATH_SEQ)

    f_error = open(PATH_DATA + "error_log.txt", "w")

    thread_pool = ThreadPool(processes=NB_THREADS)

    res = thread_pool.map(_download, urls)

    print("######## errors founds:")
    for err in res:
        if err:
            print(err)
            f_error.write('{0}\n'.format(err))

def get_urls():
    """
    get download addresses as GSM id according to the following template:
    169. TuMP2-10b
    Organism:	Mus musculus
    Source name:	mouse pancreatic tumor
    Platform: GPL15907 Series: GSE51372
    FTP download: SRA SRX364871
                            ftp://ftp-trace.ncbi.nlm.nih.gov/
                            sra/sra-instant/reads/ByExp/sra/SRX/SRX364/SRX364871/
    Sample		Accession: GSM1243834	ID: 301243834
    """
    regex_url = re.compile("(?<=!Sample_supplementary_file_1 = )ftp://.+(?=\n)")
    regex_gsm = re.compile("(?<=\^SAMPLE = )GSM[0-9]+")

    read = open(PATH_SOFT, 'r').read()

    addresses = regex_url.findall(read)
    gsms = regex_gsm.findall(read)

    return zip(gsms, addresses)

if __name__ == "__main__":
    main()
