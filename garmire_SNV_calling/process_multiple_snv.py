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

from garmire_SNV_calling.process_snv_calling import ProcessSNVCalling
from garmire_SNV_calling.config import NB_PROCESS_SNV as NB_PROCESS
from garmire_SNV_calling.config import PATH_OUTPUT
from garmire_SNV_calling.config import OUTPUT_PATH_SNV



def main():
    res = raw_input(
        "==> ready to launch SNV calling pipeline on {0} processes\n continue? (Y/n)"\
              .format(NB_PROCESS))
    if res != 'Y':
        print 'abord'
        return

    mp_analysis = Mp_Analysis()
    mp_analysis.run()

class Mp_Analysis():
    def __init__(self):
        """ """

        self.mp_queue = Queue()

        output_star = listdir(PATH_OUTPUT + "star/")

        for fil in output_star:
            if not isfile(PATH_OUTPUT + "star/" + fil + \
                          "/Aligned.sortedByCoord.out.bam"):
                print 'no star bam file for {0} skipping'.format(fil)

            if isfile("{0}/data/{1}/snv_filtered.vcf"\
                      .format(OUTPUT_PATH_SNV, fil)):
                print 'VCF file output already exists for {0} skipping...'\
                    .format(fil)
                continue

            print "file to be processed:", fil
            self.mp_queue.put(fil)

        print "\n #### now launching multiprocessing analysis #### \n"

        self.processes = [TrSNVMultiprocessing(self.mp_queue)
                          for _ in range(NB_PROCESS)]
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
    def __init__(self, input_queue):
        self.input_queue = input_queue
        self.id = randint(0, 1000000)
        self.process_snv = ProcessSNVCalling(id=self.id)
        Process.__init__(self)

    def run(self):
        while self.input_queue.qsize():
            try:
                patient = self.input_queue.get(True, 0.2)
            except Exception as e:
                print "exception:{0}".format(e)
                continue
            else:
                print "processing for file {0} with id {1}"\
                    .format(patient, self.id)
                self.process_snv.process(patient)


if __name__ == "__main__":
    main()
