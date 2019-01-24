#!/usr/bin/python

from fnmatch import fnmatch

from os import popen
from os import listdir
from os import mkdir
from os.path import isdir
from os.path import isfile
from os.path import getsize

from time import sleep
from random import random

from sys import stdout

from distutils.dir_util import mkpath

from garmire_SNV_calling.config import FASTQ_PATH
from garmire_SNV_calling.config import OUTPUT_PATH_STAR
from garmire_SNV_calling.config import PATH_STAR_SOFTWARE \
    as PATH_SOFTWARE
from garmire_SNV_calling.config import STAR_INDEX_PATH
from garmire_SNV_calling.config import STAR_INDEX_READ_LENGTH
from garmire_SNV_calling.config import STAR_THREADS as THREADS
from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN


sleep(2 * random())

if not isdir(OUTPUT_PATH_STAR):
    mkpath(OUTPUT_PATH_STAR)


def main():
    for fil in listdir(FASTQ_PATH):

        if isfile(FASTQ_PATH + fil):
            continue

        if PATTERN and not fnmatch(fil, PATTERN):
            continue

        print("====> file to be aligned:", fil)

        if not isdir(OUTPUT_PATH_STAR + fil):
            mkdir(OUTPUT_PATH_STAR + fil)

        if isfile(OUTPUT_PATH_STAR + fil + "/Aligned.sortedByCoord.out.bam") \
           and getsize(OUTPUT_PATH_STAR + fil + "/Aligned.sortedByCoord.out.bam"):
            print('bam file result alreay exists for:{0}\nskipping...'\
                .format(fil))
            continue

        fastq_str = ""

        for fastq_fil in listdir(FASTQ_PATH + fil):
            print(fastq_fil)
            if fnmatch(fastq_fil, "*.fastq"):
                fastq_str += "{0}{1}/{2} ".format(FASTQ_PATH, fil, fastq_fil)

        if not fastq_str:
            print('no fastq file found for:{0}!\nskipping'.format(fil))
            continue

        star_index_path = "{0}READ{1}".format(STAR_INDEX_PATH.rstrip('/'),
                                              STAR_INDEX_READ_LENGTH)

        cmd = "{0} --readFilesIn {1} --runThreadN {2}"\
              " --twopassMode Basic --outSAMtype BAM SortedByCoordinate" \
              "  --outFileNamePrefix {3} --genomeDir {4}"\
              .format(PATH_SOFTWARE,
                      fastq_str,
                      THREADS,
                      OUTPUT_PATH_STAR + fil + "/",
                      star_index_path
              )

        res = popen(cmd)
        c = res.read(1)

        while c:
            stdout.write(c)
            stdout.flush()
            c = res.read(1)

if __name__ == "__main__":
    main()
