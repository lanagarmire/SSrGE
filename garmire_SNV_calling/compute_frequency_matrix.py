#! /usr/bin/python

from os.path import isfile
from os.path import isdir
from os import mkdir
from os import listdir
from os import popen
from distutils.dir_util import mkpath
from sys import argv

from garmire_SNV_calling.config import FEATURE_COUNT as SOFTWARE_PATH
from garmire_SNV_calling.config import ANNOTATION_PATH
from garmire_SNV_calling.config import MATRIX_OUTPUT_PATH as OUTPUT_PATH
from garmire_SNV_calling.config import OUTPUT_PATH_STAR as STAR_PATH


############ VARIABLE ################
DEFAULT_ALIGNER = 'STAR'

if len(argv) > 1:
    DEFAULT_ALIGNER = argv[1]

OUTPUT_FILENAME = {
    'STAR': 'Aligned.sortedByCoord.out.bam',
}

PATH_DICT = {
    'STAR': STAR_PATH
}
######################################


def main():
    if DEFAULT_ALIGNER not in OUTPUT_FILENAME.keys():
        raise Exception('{0} not a regular aligner!'\
                        .format(DEFAULT_ALIGNER))

    aligner_path=PATH_DICT[DEFAULT_ALIGNER]
    output_filename = OUTPUT_FILENAME[DEFAULT_ALIGNER]

    do_expression_profile(aligner_path, output_filename)

def do_expression_profile(aligner_path, output_filename):
    """
    compute expression matrix according to aligner path results
    and output_filename (ex: output.bam)
    """

    for folder in listdir(aligner_path):

        if not isdir(aligner_path + folder):
            print 'not a folder! continuing', folder
            continue

        bam_file = "{0}/{1}/{2}"\
                   .format(aligner_path, folder, output_filename)

        if not isfile(bam_file):
            print 'no bam file for {0}'.format(bam_file)
            continue

        out_path = "{0}/{1}/{2}"\
                   .format(OUTPUT_PATH, DEFAULT_ALIGNER, folder)

        if not isdir(out_path):
            mkpath(out_path)


        cmd = "{0} -pPBCM --primary -T 1 -a {1} -o {2}/matrix_counts.txt"\
              " {3}".format(SOFTWARE_PATH,
                            ANNOTATION_PATH,
                            out_path,
                            bam_file)
        popen(cmd).read()


if __name__ == "__main__":
    main()
