#! /usr/bin/python

""" process multiple bam file with SNV"""

from multiprocessing import Process
from multiprocessing import Queue

from random import randint
from time import sleep
from os import listdir
from os.path import isfile
from os.path import isdir
from shutil import rmtree as rmdir

from garmire_SNV_calling.process_snv_GATK import ProcessGATKSNV
from garmire_SNV_calling.process_freebayes import ProcessFreebayesCaller

from garmire_SNV_calling.config import PATH_OUTPUT

from garmire_SNV_calling.config import OUTPUT_PATH_GATK
from garmire_SNV_calling.config import OUTPUT_PATH_FREEBAYES

from sys import argv


######## VARIABLE ##############################
CLEANING_MODE = True


if '--freebayes' in argv or '--do_both_callers' in argv:
    SNVCLASS = ProcessFreebayesCaller
    OUTPUT_PATH_SNV = OUTPUT_PATH_FREEBAYES
    print('GATK SNV caller used. To use Freebayes, add --freebayes option')
else:
    SNVCLASS = ProcessGATKSNV
    print('freebayes SNV caller used')
    OUTPUT_PATH_SNV = OUTPUT_PATH_GATK

if '--limit'  in argv:
    LIMIT = int(argv[argv.index('--limit') + 1 ])
else:
    LIMIT = None

if "--nb_processes" in argv:
    NB_PROCESS = eval(argv[
        argv.index("--nb_processes") + 1])
else:
    from garmire_SNV_calling.config import NB_PROCESS_SNV as NB_PROCESS
################################################


def main():
    res = raw_input(
        "==> ready to launch SNV on {0} processes\n continue? (Y/n)"\
              .format(NB_PROCESS))
    if res != 'Y':
        print('abord')
        return

    mp_analysis = Mp_Analysis()
    mp_analysis.run()


class Mp_Analysis():
    def __init__(self):
        """ """

        self.mp_queue = Queue()

        output_star = listdir(PATH_OUTPUT + "star/")

        if LIMIT:
            output_star = output_star[:LIMIT]

        for fil in output_star:
            if not isfile(PATH_OUTPUT + "star/" + fil + \
                          "/Aligned.sortedByCoord.out.bam"):
                print('no star bam file for {0} skipping'.format(fil))

                if isdir(PATH_OUTPUT + "star/" + fil) and CLEANING_MODE:
                    rmdir(PATH_OUTPUT + "star/" + fil)
                continue

            if isfile("{0}/data/{1}/snv_filtered.vcf"\
                      .format(OUTPUT_PATH_SNV, fil)):
                print('VCF file output already exists for {0} skipping...'\
                    .format(fil))
                continue

            print("file to be processed:", fil)
            self.mp_queue.put(fil)

        print("\n #### now launching multiprocessing analysis #### \n")

        self.processes = [TrSNVMultiprocessing(self.mp_queue, id=i)
                          for i in range(NB_PROCESS)]
    def _run(self):
        for p in self.processes:
            p.start()

        while self.mp_queue.qsize():
            for p in self.processes:
                if p.exitcode:
                    raise KeyboardInterrupt
            sleep(1)

    def run(self):
        try:
            self._run()

        except KeyboardInterrupt:
            for p in self.processes:
                p.terminate()


class TrSNVMultiprocessing(Process):
    """
    Launch and control several instance of SNV process
    """
    def __init__(self, input_queue, id):
        self.input_queue = input_queue
        self.id = id
        Process.__init__(self)
        self.process_snv = SNVCLASS(id=self.id)

    def run(self):
        while self.input_queue.qsize():
            try:
                patient = self.input_queue.get(True, 0.2)
            except Exception as e:
                print("exception:{0}".format(e))
                continue
            else:
                print("processing for file {0} with id {1}"\
                    .format(patient, self.id))

                if '--do_both_callers' in argv:
                    error = self.process_snv.process_ALL_callers(patient)
                else:
                    error = self.process_snv.process(patient)

                if error is not None:
                    print('error {1} found for patient: {0}'.format(patient, error))


if __name__ == "__main__":
    main()
