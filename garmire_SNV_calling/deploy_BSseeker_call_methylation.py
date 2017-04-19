""" """

from multiprocessing import Pool

from os import popen

from os.path import isdir
from os.path import isfile
from os.path import split as pathsplit

from glob import glob

from time import sleep
from random import random

from fnmatch import fnmatch

from distutils.dir_util import mkpath

from garmire_SNV_calling.config import PATH_OUTPUT
from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN
from garmire_SNV_calling.config import BSSEEKER2_REP
from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import PYTHON


################ VARIABLE ##################################

BAM_PATH = PATH_OUTPUT + '/BSseeker/'
PROCESS_THREADS = 2
BISMARK_OPTION = ''
REF_GENOME_PATH = pathsplit(REF_GENOME)[0]
REF_GENOME_RRBS_DB = REF_GENOME_PATH + '/genome.fa_rrbs_20_500_bowtie2/'
############################################################


sleep(2 * random())
if not isdir(BAM_PATH):
    mkpath(BAM_PATH)


def main():
    pool = Pool(PROCESS_THREADS)
    # process_one_file(glob(BAM_PATH + '/*')[0])
    pool.map(process_one_file, glob(BAM_PATH + '/*'))

def process_one_file(folder):
    """ """
    print(folder)
    if isfile(folder):
        return False

    if PATTERN and not fnmatch(folder, PATTERN):
        return False

    print("====> folder to be processed:", folder)

    input_bam_file_name = glob(folder + '/*.bam')

    if not input_bam_file_name:
        print('no bam file detected for :{0}\nskipping...'\
            .format(folder))
        return False

    if len(input_bam_file_name) > 1:
        print('multiple bam files detected: {0}. selecting the first'.format(
            input_bam_file_name))

    input_bam = input_bam_file_name[0]
    output_file = input_bam.rsplit('.', 1)[0] + '.CpG.CGmap'

    cmd = "{0} {1}/bs_seeker2-call_methylation.py -i {2}  --CGmap {3} --db {4} --txt " \
        .format(PYTHON,
                BSSEEKER2_REP,
                input_bam,
                output_file,
                REF_GENOME_RRBS_DB,
        )

    _run_cmd(cmd)
    _run_cmd('rm {0}/*_sorted*'.format(folder))

    return True


def _run_cmd(cmd, *args):
    """run cmd"""
    popen(cmd).read()


if __name__ == "__main__":
    main()
