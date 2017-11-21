""" """

from multiprocessing import Pool

from os import popen

from os import listdir
from os import mkdir
from os.path import isdir
from os.path import isfile
from os.path import getsize
from os.path import split as pathsplit

import re

from glob import glob

from time import sleep
from random import random

from fnmatch import fnmatch

from distutils.dir_util import mkpath

from garmire_SNV_calling.config import FASTQ_PATH
from garmire_SNV_calling.config import PATH_OUTPUT
from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN
from garmire_SNV_calling.config import BSSEEKER2_REP
from garmire_SNV_calling.config import PYTHON
from garmire_SNV_calling.config import BOWTIE_REP
from garmire_SNV_calling.config import DO_TRIMGALORE
from garmire_SNV_calling.config import TRIMGALORE_REP


################ VARIABLE ##################################

OUTPUT_PATH = PATH_OUTPUT + '/BSseeker/'
PROCESS_THREADS = 2
BISMARK_OPTION = ''
REF_GENOME_PATH = pathsplit(REF_GENOME)[0]
############################################################


sleep(2 * random())
if not isdir(OUTPUT_PATH):
    mkpath(OUTPUT_PATH)


def main():
    pool = Pool(PROCESS_THREADS)
    pool.map(process_one_file, listdir(FASTQ_PATH))

def process_one_file(fil):
    """ """
    print(fil)
    if isfile(FASTQ_PATH + fil):
        return False

    if PATTERN and not fnmatch(fil, PATTERN):
        return False

    print("====> file to be aligned:", fil)

    if not isdir(OUTPUT_PATH + fil):
        mkdir(OUTPUT_PATH + fil)

    bam_file_name = glob(OUTPUT_PATH + fil + '/*.bam')

    if bam_file_name \
       and getsize(bam_file_name[0]):
        print('bam file result alreay exists for:{0}\nskipping...'\
            .format(bam_file_name[0]))
        return False

    fastq_str = ""

    fastq_files = list(set(glob(FASTQ_PATH + fil + '/*.fastq')))
    print('fastq files founds: {0}'.format(fastq_files))

    stdout = open(OUTPUT_PATH + fil + "/log.out", 'w')

    if len(fastq_files) > 2:
        print('tow many fastq files!')
        return False

    elif len(fastq_files) == 2:
        fastq_1 = [fastq for fastq in fastq_files
                   if re.match('.+_1\.fastq', fastq)]

        assert(fastq_1)

        fastq_1 = fastq_1[0]

        fastq_2 = [fastq for fastq in fastq_files
                   if re.match('.+_2\.fastq', fastq,)]
        assert(fastq_2)

        fastq_2 = fastq_2[0]

        if DO_TRIMGALORE:
            cmd_trim = "{0}/trim_galore {1} {2} --paired --no_report_file -o {3}".format(
                TRIMGALORE_REP, fastq_1, fastq_2, FASTQ_PATH + fil)
            _run_cmd(cmd_trim, stdout)

            fastq_1 = '{0}_val_1.fq'.format(fastq_1.rsplit('.', 1)[0])
            fastq_2 = '{0}_val_2.fq'.format(fastq_2.rsplit('.', 1)[0])

        fastq_str = ' -1 {0} -2 {1} '.format(fastq_1, fastq_2)

    elif len(fastq_files) == 1:
        if DO_TRIMGALORE:
            fastq_file = fastq_files[0]
            cmd_trim = "{0}/trim_galore {1} --no_report_file -o {2}".format(
                TRIMGALORE_REP, fastq_file, FASTQ_PATH + fil)
            _run_cmd(cmd_trim, stdout)
            fastq_file = '{0}_trimmed.fq'.format(fastq_file.rsplit('.', 1)[0])

        fastq_str = ' -i {0} '.format(fastq_file)

    if not fastq_str:
        print('no fastq file found for:{0}!\nskipping'.format(fil))
        return False

    cmd = "{0} {1}/bs_seeker2-align.py -g {2}"\
          " --aligner=bowtie2 -p {3} --db {4} -r {5}"\
        .format(PYTHON,
                BSSEEKER2_REP,
                REF_GENOME,
                BOWTIE_REP,
                REF_GENOME_PATH,
                fastq_str
    )

    _run_cmd(cmd, stdout)
    _run_cmd('mv {0}/*.bam {1}/; mv {0}/*log {1}/ '.format(
        FASTQ_PATH + fil, OUTPUT_PATH + fil), stdout)

    return True


def _run_cmd(cmd, *args):
    """run cmd"""
    popen(cmd).read()


if __name__ == "__main__":
    main()
