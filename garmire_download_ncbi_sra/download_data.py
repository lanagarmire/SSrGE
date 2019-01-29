#! /usr/bin/python

"""
download data from NCBI according to GEO accession number
"""

from os.path import isdir
from os.path import isfile
from os import popen

from os import mkdir

import urllib2

from distutils.dir_util import mkpath

import json

from garmire_download_ncbi_sra.config import PATH_DATA
from garmire_download_ncbi_sra.config import NB_THREADS
from garmire_download_ncbi_sra.config import LIMIT

from garmire_SNV_calling.bash_utils import exec_cmd

from urllib2 import URLError

import re

from multiprocessing.pool import ThreadPool

from time import sleep


############ VARIABLES ############
PATH_SEQ = PATH_DATA + '/fastq/'

if not isdir(PATH_SEQ):
    mkdir(PATH_SEQ)
###################################


def main():
    download_data()

def _download_old(url):
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

def _download(data, verbose=True):
    """ """
    gsm, url_address = data

    waiting_list = [10, 20, 30]

    f_name = "{0}/{1}.sra".format(PATH_SEQ, gsm)

    if isfile('{0}/download_successfull.log'.format(PATH_SEQ)):
        msg = 'file {0} already downloaded. skipping...'.format(f_name)
        print(msg)
        return msg
    while True:
        try:
            url = urllib2.urlopen(url_address).read()
        except URLError as e:
            if waiting_list:
                sleep_time = waiting_list.pop()
                print('error when downloading: {1} sleeping {0} s...'.format(sleep_time, e))
                sleep(sleep_time)
            else:
                raise e
        else:
            break

    srr = re.findall('run=(?P<srr>SRR[0-9]+)', url)[0]

    srr_url = "ftp://ftp-trace.ncbi.nlm.nih.gov"\
              "/sra/sra-instant/reads/ByRun/sra/SRR/{0}/{1}/{1}.sra".format(
                  srr[0:6], srr)

    print('downloading: {0}'.format(srr_url))

    if verbose:
        verb = ''
    else:
        verb = '--no-verbose'
    cmd = "wget {2} -O {0} {1} ".format(f_name, srr_url, verb)
    print('launching: {0}'.format(cmd))

    exec_cmd(cmd)

    msg = '{0} successfully downloaded'.format(f_name)
    print(msg)

    f_log = open('{0}/download_successfull.log'.format(PATH_SEQ), 'w')
    f_log.write('download complete')

    return msg

def download_data():
    """download dataset from ncbi """

    urls = get_urls()

    if LIMIT:
        urls = urls[:LIMIT]

    if not isdir(PATH_SEQ):
        mkpath(PATH_SEQ)

    f_error = open(PATH_DATA + "/error_log.txt", "w")

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
    f_meta = open('{0}/metadata.json'.format(PATH_DATA))
    metadata = json.load(f_meta)

    gsms, urls = [], []

    for sample in metadata:
        if 'SRA' in metadata[sample]:
            gsms.append(sample)
            urls.append(metadata[sample]['SRA'])

    return zip(gsms, urls)

if __name__ == "__main__":
    main()
