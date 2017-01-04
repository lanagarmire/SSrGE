#! /usr/bin/python

"""load data """

from garmire_SSrGE.config import PROJECT_PATH
from garmire_SSrGE.config import SOFT_PATH

from garmire_SSrGE.generate_refgenome_index import INDEX_SAVE_PATH
from garmire_SSrGE.generate_refgenome_index import main as generate_refgenome

import cPickle
import re

from time import time
from os.path import isfile

from collections import defaultdict


def load_indexes(path_indexes=INDEX_SAVE_PATH):
    t = time()

    if not isfile(path_indexes + 'index_start.pickle'):
        print 'indexes not found. Creating indexes...'
        generate_refgenome()

    with open(path_indexes + 'index_start.pickle', 'r') as f:
        index_start = cPickle.load(f)
    with open(path_indexes + 'index_end.pickle', 'r') as f:
        index_end = cPickle.load(f)
    with open(path_indexes + 'position_index.pickle', 'r') as f:
        position_index = cPickle.load(f)

    print 'gene position indexes loaded in {0} s'.format(time() - t)
    return index_start, index_end, position_index

def process_line_from_vcf_file(line):
    """ process one line from the svf file"""

    if line[0] == '#':
        return
    line = line.split('\t')

    snv_id = None

    # process only passed SNV
    if line[6] != 'PASS':
        return

    # take annotation
    if line[2] != '.':
        snv_id = line[2]

    chrid, start = line[0], int(line[1])
    end = start
    return chrid, start, end, snv_id

def load_gsm_and_sample_names_from_soft(soft_path=SOFT_PATH):
    """
    load GSM and sample names from soft

    return:
        dict(GSM:sample name)
    """
    if not soft_path:
        return defaultdict(str)

    regex_gsm = re.compile("(?<=\^SAMPLE = )GSM[0-9]+")
    regex_name = re.compile("(?<=!Sample_title = ).+(?=\n)")

    if not isfile(soft_path):
        return {}

    read = open(soft_path, 'r').read()

    gsms = regex_gsm.findall(read)
    names = regex_name.findall(read)

    return defaultdict(str, zip(gsms, names))

def process_line_from_annotated_vcf_file(line):
    """ process one line from the vcf file"""

    if line[0] == '#':
        return
    line = line.split('\t')

    # process only passed SNV
    if line[6] != 'PASS':
        return

    if ONLY_KNOWN_SNV:
        # process only SNV with a known id
        if line[2] == '.':
            return
    results = []
    annotations = line[7].split(';')[-1].split(',')

    for annotation in annotations:
        annotation = annotation.split('|')

        if len(annotation) < 7:
            continue
        if not annotation[1]:
            continue

        result = {'type': annotation[1],
                  'impact': annotation[2],
                  'gene impacted': annotation[3],
                  'feature type': annotation[5],
                  'biotype': annotation[7]
        }
        results.append(result)
    chrid, start = line[0], int(line[1])
    return (annotation[3], start), results
