from config import PATH_OUTPUT
from config import OUTPUT_PATH_STAR
from config import MONOVAR_REP
from config import MONOVAR_SAMTOOLS
from config import REF_GENOME
from config import PYTHON
from config import NB_PROCESS_SNV

from os.path import isdir
from os import mkdir

from glob import glob

from random import sample

from os import remove
from os import popen

from multiprocessing import Pool


######## LOCAL VARIABLES ############################
PATH_MONOVAR = '{0}/monovar/'.format(PATH_OUTPUT)
CHUNCK_SIZE = 20
THREAD_NB = 3
#####################################################


if not isdir(PATH_MONOVAR):
    mkdir(PATH_MONOVAR)


def main():
    """ """
    create_list_file()
    launch_monovar()

def launch_monovar():
    """
    """
    cmd_list = []

    for fil in glob('{0}/monovar_input*.txt'.format(PATH_MONOVAR)):

        cmd = '{0} mpileup -BQ0 -d10000 -f {1} -b {6} '\
              '| {3} {4}/src/monovar.py -p 0.002 -a 0.2 -t 0.05 -f {1}'\
              ' -b {6} -m {5} -o {7}.vcf'\
              .format(MONOVAR_SAMTOOLS,
                      REF_GENOME,
                      PATH_MONOVAR,
                      PYTHON,
                      MONOVAR_REP,
                      NB_PROCESS_SNV,
                      fil, fil.rsplit('.', 1)[0])

        cmd_list.append(cmd)
    # import ipdb;ipdb.set_trace()

    pool = Pool(THREAD_NB)
    pool.map(_multiprocessing_func, cmd_list)

def _multiprocessing_func(cmd):
    """ """
    print('###### command launched:\n{0}\n########'.format(cmd))
    popen(cmd).read()

def create_list_file():
    """
    create the input file used by monovar containing all the input files
    """
    for fil in glob('{0}/monovar_input*'.format(PATH_MONOVAR)):
        remove(fil)

    file_list = set()

    for folder in glob('{0}/*'.format(OUTPUT_PATH_STAR)):
        if not isdir(folder):
            continue

        file_list.add('{0}/Aligned.sortedByCoord.out.bam'.format(folder))

    nb_file = 0

    chunck_size = CHUNCK_SIZE

    while file_list:
        if len(file_list) < chunck_size:
            chunck_size = len(file_list)

        sample_list = sample(file_list, chunck_size)
        file_list = file_list.difference(sample_list)

        f_input = open('{0}/monovar_input_{1}.txt'.format(PATH_MONOVAR, nb_file), 'w')

        for fil in sample_list:
            f_input.write('{0}\n'.format(fil))

            nb_file += 1


if __name__ == "__main__":
    main()
