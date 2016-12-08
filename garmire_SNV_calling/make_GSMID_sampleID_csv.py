#!/usr/bin/python

"""
read soft summary file from GEO page and map GSE to cell name
"""

from garmire_SNV_calling.config import GLOBAL_DATA_ROOT
from garmire_SNV_calling.config import PROJECT_NAME
from garmire_SNV_calling.config import SOFT_PATH

from os.path import isfile

import re

def main():
    csv_path = "{0}/{1}/{1}.csv"\
               .format(GLOBAL_DATA_ROOT, PROJECT_NAME)

    if not isfile(SOFT_PATH):
        print "error! no file: {0}".format(SOFT_PATH)
        return 1

    f_soft = open(SOFT_PATH, 'r').read()
    f_csv = open(csv_path, 'w')

    gse_list = re.findall("(?<=\^SAMPLE \= )\w+", f_soft)
    id_list = re.findall("(?<=!Sample_title \= ).+(?!\n)", f_soft)

    for gse, ids in zip(gse_list, id_list):
        ids = ids.replace(' ', '_')
        f_csv.write("{0};{1}\n".format(gse, ids))

    print "done"

if __name__ == "__main__":
    main()
