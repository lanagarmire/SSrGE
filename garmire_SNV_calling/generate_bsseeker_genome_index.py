"""generate STAR GENOME INDEX"""

from sys import stdout as sys_stdout
from os import popen
from os import mkdir
from os.path import isdir
from os.path import split as pathsplit

from distutils.dir_util import mkpath

from garmire_SNV_calling.config import BSSEEKER2_REP
from garmire_SNV_calling.config import BSSEQ_INDEX_PATH
from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import PYTHON
from garmire_SNV_calling.config import BOWTIE_REP


################ VARIABLE ################
REF_GENOME_PATH = pathsplit(REF_GENOME)[0]
##########################################


def main():
    """ """
    bsseq_index_path = BSSEQ_INDEX_PATH
    print "######## computing BS-seq index ########\npath:{0}\n"\
        .format(bsseq_index_path)

    if not isdir(bsseq_index_path):
        mkpath(bsseq_index_path)

    cmd = "{0} {1}/bs_seeker2-build.py -f {2}"\
          " --aligner=bowtie2 -p {3} --db {4} -r"\
        .format(PYTHON,
                BSSEEKER2_REP,
                REF_GENOME,
                BOWTIE_REP,
                REF_GENOME_PATH
    )

    stdout = popen(cmd)
    c = stdout.read(1)

    while c:
        sys_stdout.write(c)
        sys_stdout.flush()
        c = stdout.read(1)

if __name__ == "__main__":
    main()
