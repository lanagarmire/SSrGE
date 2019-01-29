#!/usr/bin/python
from fnmatch import fnmatch

from os import listdir
from os import mkdir
from os.path import isdir
from os.path import isfile
from os.path import getsize

from garmire_SNV_calling.bash_utils import exec_cmd

from time import sleep
from random import random
from sys import argv

from os import popen

from distutils.dir_util import mkpath

from garmire_SNV_calling.config import PATH_STAR_SOFTWARE \
    as PATH_SOFTWARE

from garmire_SNV_calling.config import STAR_THREADS as THREADS

from garmire_SNV_calling.bash_utils import printf


############ VARIABLES ############################################

from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN
from garmire_SNV_calling.config import FASTQ_PATH
from garmire_SNV_calling.config import STAR_INDEX_PATH
from garmire_SNV_calling.config import STAR_INDEX_READ_LENGTH

from garmire_SNV_calling.config import OUTPUT_PATH_STAR \
    as OUTPUT_PATH

from garmire_SNV_calling.config import SIMULATED_REF_GENOME

###################################################################


def star_analysis(
        output_path=OUTPUT_PATH,
        fastq_path=FASTQ_PATH,
        pattern=PATTERN,
        star_index_path=STAR_INDEX_PATH,
        star_index_read_length=STAR_INDEX_READ_LENGTH,
        simulated_ref_genome=SIMULATED_REF_GENOME,
        path_software=PATH_SOFTWARE,
        threads=THREADS,
        cufflinks_compatibility=None,
        custom_star_index_name=True,
        stdout=None,
        printf=printf):
    """
    """
    sleep(2 * random())

    options = ''

    if cufflinks_compatibility:
        options = '--outSAMstrandField intronMotif'\
                  ' --outFilterIntronMotifs RemoveNoncanonical'

    if not isdir(output_path):
        mkpath(output_path)

    for fil in listdir(fastq_path):
        if isfile(fastq_path + '/' + fil):
            continue

        if pattern and not fnmatch(fil, pattern):
            continue

        printf("====> file to be aligned: {0}".format(fil))

        if not isdir(output_path + fil):
            mkdir(output_path + fil)

        if isfile(output_path + fil + "/Aligned.sortedByCoord.out.bam") \
           and getsize(output_path + fil + "/Aligned.sortedByCoord.out.bam"):
            printf('bam file result alreay exists for:{0}\nskipping...'\
                .format(fil))
            continue

        fastq_str = ""

        for fastq_fil in sorted(listdir(fastq_path + '/' +  fil)):
            if fnmatch(fastq_fil, "*.fastq"):
                fastq_str += "{0}/{1}/{2} ".format(fastq_path, fil, fastq_fil)

        if not fastq_str:
            printf('no fastq file found for:{0}!\nskipping'.format(fil))
            continue

        if custom_star_index_name:
            star_index_path_ready = "{0}READ{1}".format(star_index_path.rstrip('/'),
                                                        star_index_read_length)

            if simulated_ref_genome:
                    star_index_path_ready = star_index_path_ready.rstrip('/') \
                          + 'SIM{0}/'.format(simulated_ref_genome)
        else:
            star_index_path_ready = star_index_path

        tmp_path = '{0}/_STARtmp'.format(output_path + '/' + fil + "/")

        if isdir(tmp_path):
            exec_cmd('rm -r {0}'.format(tmp_path), stdout)

        cmd = "{0} --readFilesIn {1} --runThreadN {2}"\
              " --twopassMode Basic --outSAMtype BAM SortedByCoordinate" \
              "  --outTmpDir {5} --outFileNamePrefix {3} --genomeDir {4} {6}"\
              .format(path_software,
                      fastq_str,
                      threads,
                      output_path + '/' + fil + "/",
                      star_index_path_ready,
                      tmp_path,
                      options
              )

        printf('star cmd to be launched:{0}'.format(cmd))
        exec_cmd(cmd, stdout)

def check_star_folder(new_star_path):
    """
    """
    if isfile('{0}/star_aligned_successfull.log'.format(new_star_path)):
        return '#### STAR already aligned successfully in: {0}'.format(new_star_path)

def clean_star_folder(path_star_results):
    """
    """
    path_bam = '{0}/Aligned.sortedByCoord.out.bam'.format(path_star_results)
    path_log_final = '{0}/Log.final.out'.format(path_star_results)
    path_sj = '{0}/SJ.out.tab'.format(path_star_results)

    assert(isfile(path_bam) and isfile(path_log_final) and isfile(path_sj))

    popen('rm -r {0}/_STAR*'.format(path_star_results)).read()
    popen('rm -r {0}/Log.out'.format(path_star_results)).read()
    popen('rm -r {0}/Log.progress.out'.format(path_star_results)).read()

    f_log = open('{0}/star_aligned_successfull.log'.format(path_star_results), 'w')
    f_log.write('STAR successfull')

    return


if __name__ == "__main__":
    star_analysis()
