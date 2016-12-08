
""" process one fastqc report"""

from os import popen
from os import listdir

from os.path import isdir
from os.path import isfile
from os.path import getsize
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

from garmire_SNV_calling.config import FASTQC
from garmire_SNV_calling.config import PATH_OUTPUT
from garmire_SNV_calling.config import FASTQ_PATH

from garmire_SNV_calling.process_multiple_generic import MPI

############ VARIABLES ############################################
SRR_TO_PROCESS = "" # for debug purpose
PROCESS_ID = randint(0, 1000000)

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
        input_list = listdir(FASTQ_PATH)
        mpi = MPI(input_list=input_list,
                  ProcessClass=ProcessFastqC,
                  nb_processes=NB_THREADS)
        mpi.run()
    else:
        process_fastqc = ProcessFastqC(id=PROCESS_ID)
        process_fastqc.process(SRR_TO_PROCESS)

class ProcessFastqC():
    """ Process Fastqc report"""
    def __init__(self,
                 path_to_data=PATH_OUTPUT,
                 id="1",
                 clean_tmp=True,
    ):
        self.output_path = PATH_OUTPUT + '/fastqc/'
        self.path_to_data = path_to_data
        self.time_start = None
        self.id = str(id)
        self.stdout = None

    def process(self, srr_to_process=SRR_TO_PROCESS):
        """
        process one fastq file using fastqc
        """
        tmppath = self.output_path + "/tmp/" + self.id
        outpath = self.output_path + "/data/"

        if not isdir(FASTQ_PATH + srr_to_process):
            print '{0} is not a folder!'.format(
                FASTQ_PATH + srr_to_process)
            return

        if isdir("{1}/{0}_fastqc"\
                 .format(srr_to_process, outpath)):
            print '{0} output already exists'.format(
                "{1}/{0}_fastqc"\
                 .format(srr_to_process, outpath))
            return

        sleep(random())
        if not isdir(tmppath):
            mkpath(tmppath)
        if not isdir(outpath):
            mkpath(outpath)

        popen("rm {0}/*".format(tmppath)).read()
        path_fastq = ""

        for fil in listdir(FASTQ_PATH + srr_to_process):
            if fnmatch(fil, '*.fastq'):
                path_fastq = '{0}/{1}/{2}'.format(FASTQ_PATH,
                                                  srr_to_process,
                                                  fil)
                fil = fil.rsplit('.', 1)[0]
                break
        if not path_fastq:
            print 'No fastq file for :{0}'.format(path_fastq)
            return

        self.stdout = open(tmppath + '/stdout.log', 'w')

        cmd = "{0} {1} -o {2} -d {2} --extract"\
              .format(FASTQC, path_fastq, tmppath)

        self._exec_cmd(cmd)
        self._exec_cmd("mv {0}/{1}_fastqc {3}/{2}_fastqc"\
                       .format(tmppath, fil, srr_to_process, outpath))
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
