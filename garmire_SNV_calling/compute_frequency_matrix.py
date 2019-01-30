#! /usr/bin/python

from os.path import isfile
from os.path import isdir

from os import listdir

from os.path import getsize
from distutils.dir_util import mkpath

from garmire_SNV_calling.config import FEATURE_COUNT
from garmire_SNV_calling.config import ANNOTATION_PATH
from garmire_SNV_calling.config import MATRIX_OUTPUT_PATH as OUTPUT_PATH
from garmire_SNV_calling.config import OUTPUT_PATH_STAR as STAR_PATH
from garmire_SNV_calling.config import STAR_THREADS

from multiprocessing import Pool

from garmire_SNV_calling.bash_utils import exec_cmd


############ VARIABLE ################
DEFAULT_ALIGNER = 'STAR'

OUTPUT_FILENAME = {
    'STAR': 'Aligned.sortedByCoord.out.bam'
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
    cmd_list = []

    for folder in listdir(aligner_path):

        if not isdir(aligner_path + folder):
            print('not a folder! continuing', folder)
            continue

        bam_file = "{0}/{1}/{2}"\
                   .format(aligner_path, folder, output_filename)

        if not isfile(bam_file):
            print('no bam file for {0}'.format(bam_file))
            continue

        out_folder = "{0}/{1}/{2}"\
                   .format(OUTPUT_PATH, DEFAULT_ALIGNER, folder)

        out_file = '{0}/{1}'.format(out_folder, "matrix_count.txt")

        if isfile(out_file) and getsize(out_file):
            print('expression matrix already exists for: {0}'.format(out_folder))
            continue

        cmd_list.append((bam_file, out_folder))

    pool = Pool(STAR_THREADS)
    pool.map(_multiprocess_func, cmd_list)


def _multiprocess_func(inp):
    """ """
    bam_file, out_folder = inp
    bam_file_to_expression_matrix(bam_file, out_folder)


def bam_file_to_expression_matrix(
        bam_file,
        out_folder,
        feature_count=FEATURE_COUNT,
        annotation_path=ANNOTATION_PATH,
        stdout=None,
        matrix_name="matrix_counts.txt"):
    """ """
    if not isdir(out_folder):
            mkpath(out_folder)

    cmd = "{0} -pPBCM --primary -T 1 -a {1} -o {2}/{4}"\
          " {3}".format(feature_count,
                        annotation_path,
                        out_folder,
                        bam_file,
                        matrix_name)
    print('launching cmd: {0}\n'.format(cmd))
    try:
        exec_cmd(cmd, stdout)
    except Exception as e:
        print('exception with featureCount cmd: {0}\n'.format(cmd))
        print('exception: {0}\n'.format(e))

        assert(isfile('{0}/{1}'.format(out_folder, matrix_name)))


if __name__ == "__main__":
    main()
