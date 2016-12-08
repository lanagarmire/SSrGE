#! /usr/bin/python

""" check overall statistics for all log files from fastqc report"""

from os import listdir
from os.path import isfile
import re

from collections import Counter

from garmire_SNV_calling.config import PATH_OUTPUT

PATH_OUTPUT_FASTQC = PATH_OUTPUT + '/fastqc/data/'


def main():
    make_aligner_quality_csv()

def make_aligner_quality_csv():
    """ """
    regex_status = "(?<=Sequence Duplication Levels\t)\w+"
    regex_status = re.compile(regex_status)

    stats_status = {}

    for folder in listdir(PATH_OUTPUT_FASTQC):
        log_file = "{0}/{1}/fastqc_data.txt"\
                   .format(PATH_OUTPUT_FASTQC, folder)

        if not isfile(log_file):
            continue

        read = open(log_file, 'r').read()
        status = regex_status.findall(read)[0]
        sample = folder.rsplit('_fastqc', 1)[0]
        stats_status[sample] = status

    f_csv = open(PATH_OUTPUT + '/deduplicated_check.csv', 'w')

    for key in stats_status:
        f_csv.write('{0};{1}\n'.format(key, stats_status[key]))


if __name__ == "__main__":
    main()
