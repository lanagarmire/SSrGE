from multiprocessing import Pool

from os import popen
from subprocess import Popen
from subprocess import PIPE

from os import listdir
from os import mkdir
from os.path import isdir
from os.path import isfile
from os.path import getsize
from os.path import split as pathsplit

import re

from sys import stdout as STDOUT

from glob import glob

from time import sleep
from random import random
from sys import argv
from sys import stdout

from fnmatch import fnmatch

from distutils.dir_util import mkpath

from garmire_SNV_calling.config import FASTQ_PATH
from garmire_SNV_calling.config import PATH_OUTPUT
from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN


################ VARIABLE ##################################
BISMARK_SOFTWARE = '/home/opoirion/prog/Bismark/bismark'
OUTPUT_PATH = PATH_OUTPUT + '/bismark/'
THREADS = 4
PROCESS_THREADS = 4
BISMARK_OPTION = ''
REF_GENOME_DIR = pathsplit(REF_GENOME)[0]
############################################################


sleep(2 * random())
if not isdir(OUTPUT_PATH):
    mkpath(OUTPUT_PATH)


def main():
    pool = Pool(PROCESS_THREADS)
    # process_one_file(listdir(FASTQ_PATH)[0])
    pool.map(process_one_file, listdir(FASTQ_PATH))

def process_one_file(fil):
    """ """
    if isfile(FASTQ_PATH + fil):
        return

    if PATTERN and not fnmatch(fil, PATTERN):
        return

    print "====> file to be aligned:", fil

    if not isdir(OUTPUT_PATH + fil):
        mkdir(OUTPUT_PATH + fil)

    bam_file_name = glob(OUTPUT_PATH + fil + '/*.bam')

    if bam_file_name \
       and getsize(bam_file_name[0]):
        print 'bam file result alreay exists for:{0}\nskipping...'\
            .format(bam_file_name[0])
        return

    fastq_str = ""

    fastq_files = list(set(glob(FASTQ_PATH + fil + '/*.fastq')))
    print 'fastq files founds: {0}'.format(fastq_files)

    if len(fastq_files) > 2:
        print 'tow many fastq files!'
        return

    elif len(fastq_files) == 2:
        fastq_1 = [fastq for fastq in fastq_files
                   if re.match('.+_1\.fastq', fastq)]

        assert(fastq_1)

        fastq_1 = fastq_1[0]

        fastq_2 = [fastq for fastq in fastq_files
                   if re.match('.+_2\.fastq', fastq,)]
        assert(fastq_2)

        fastq_2 = fastq_2[0]
        fastq_str = ' -1 {0} -2 {1} '.format(fastq_1, fastq_2)

    elif len(fastq_files) == 1:
        fastq_str = ' {0} '.format(fastq_files[0])

    stdout = open(OUTPUT_PATH + fil + "/log.out", 'w')

    if not fastq_str:
        print 'no fastq file found for:{0}!\nskipping'.format(fil)
        return

    cmd = "{0} -p {1} -o {2} --temp_dir {2} {3} --genome {4} {5} > {2}/stdlog.out" \
          .format(BISMARK_SOFTWARE,
                  THREADS,
                  OUTPUT_PATH + fil + "/",
                  BISMARK_OPTION,
                  REF_GENOME_DIR,
                  fastq_str
          )

    _run_cmd(cmd, stdout)


def _run_cmd(cmd, stdout):
    """run cmd"""

    process = popen(cmd).read()


if __name__ == "__main__":
    main()
