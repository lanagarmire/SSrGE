
""" process one fastqc report"""

from os import popen
from os import listdir

from os.path import isdir
from os.path import isfile
from subprocess import Popen

from distutils.dir_util import mkpath
from shutil import rmtree
from shutil import copyfile
from shutil import move
from sys import stdout as STDOUT
from sys import argv
from random import randint
from random import random
from time import sleep
from time import time
from fnmatch import fnmatch

from garmire_SNV_calling.config import OUTPUT_PATH_SNV
from garmire_SNV_calling.config import SNPEFF
from garmire_SNV_calling.config import JAVA
from garmire_SNV_calling.config import SNPEFF_DB

from garmire_SNV_calling.process_multiple_generic import MPI

############ VARIABLES ############################################
SRR_TO_PROCESS = "" # for debug purpose
PROCESS_ID = randint(0, 1000000)
INPUT_PATH = OUTPUT_PATH_SNV + '/data/'

if "--specific_folder" in argv:
    SRR_TO_PROCESS = argv[
        argv.index("--specific_folder") + 1]
if "--process_id" in argv:
    PROCESS_ID = int(argv[
        argv.index("--process_id") + 1])
if "--nb_threads" in argv:
    NB_THREADS = int(argv[
        argv.index("--nb_threads") + 1])
else:
    NB_THREADS = None
###################################################################


def main():
    if NB_THREADS:
        input_list = listdir(INPUT_PATH)
        mpi = MPI(input_list=input_list,
                  ProcessClass=ProcessAnnotateSNV,
                  nb_processes=NB_THREADS)
        mpi.run()
    else:
        process_annotate_snv = ProcessAnnotateSNV(id=PROCESS_ID)
        process_annotate_snv.process(SRR_TO_PROCESS)

class ProcessAnnotateSNV():
    """
    Process SNV annotation using snpEff software
    """
    def __init__(self,
                 path_to_data=OUTPUT_PATH_SNV,
                 id="1",
                 clean_tmp=True,
    ):
        self.path_to_data = path_to_data
        self.time_start = None
        self.id = str(id)
        self.stdout = None

    def process(self, srr_to_process=SRR_TO_PROCESS):
        """
        process one fastq file using fastqc
        """
        tmppath = self.path_to_data + "/tmp/" + self.id
        inputpath = self.path_to_data + "/data/"
        input_file = '{0}/{1}/snv_filtered.vcf'\
                     .format(inputpath, srr_to_process)
        tmp_file = '{0}/snv_filtered_annotated.vcf'\
                      .format(tmppath)
        output_file = '{0}/{1}/snv_filtered_annotated.vcf'\
                      .format(inputpath, srr_to_process)

        if not isdir(inputpath):
            print '{0} is not a folder!'.format(
                self.path_to_data + srr_to_process)
            return

        if not isfile(input_file):
            print 'no input file: {0}!'.format(
                input_file)
            return

        if isfile(output_file):
            print '{0} already exists!'.format(
                output_file)
            return

        sleep(random())
        if not isdir(tmppath):
            mkpath(tmppath)

        popen("rm {0}/*".format(tmppath)).read()

        self.stdout = open(tmppath + '/stdout.log', 'w')

        cmd = "{0} -jar {1} eff -v {2} {3} -noStats > {4}"\
              .format(JAVA, SNPEFF, SNPEFF_DB, input_file, tmp_file)

        self._exec_cmd(cmd)
        self._exec_cmd("mv {0} {1}"\
                       .format(tmp_file, output_file))
        self._exec_cmd("mv {0}/stdout.log {1}/{2}/snv_annotation.log"\
                       .format(tmppath, inputpath, srr_to_process))
        self.stdout.close()
        popen('rm -r {0}/*'.format(tmppath))

    def _exec_cmd(self, cmd):
        """ execute cmd """
        process = Popen(cmd,
                        stdout=self.stdout,
                        stderr=self.stdout,
                        shell=True)

        process.communicate()
        if process.returncode:
            raise Exception('{0} raise non 0 return code!\n'\
                            .format(cmd))


if __name__ == "__main__":
    main()
