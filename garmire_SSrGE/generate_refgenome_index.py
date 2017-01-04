#! /usr/bin/python

""" Generate homemade reference genome gtf index using python objects """

from collections import defaultdict
import re
from time import time
import cPickle

from os.path import isdir
from distutils.dir_util import mkpath

from garmire_SSrGE.config import INDEX_SAVE_PATH
from garmire_SSrGE.config import GTF_PATH


def main():
    t = time()
    print 'loading index...'
    index_start, index_end = load_indexed_gene_annotations()
    position_index = create_position_indexes(index_start, index_end)
    print 'done in {0} s'.format(time() - t)
    save_indexes(index_start, index_end, position_index)

    return True

def save_indexes(index_start,
                 index_end,
                 position_index,
                 save_path=INDEX_SAVE_PATH):
    """ """
    if not isdir(save_path):
        r = raw_input("{0} doesn't exist create it? (y/N)".format(save_path))
        if r != 'y':
            return

    mkpath(save_path)

    with open(save_path + 'index_start.pickle', 'w') as f:
        cPickle.dump(index_start, f)
    with open(save_path + 'index_end.pickle', 'w') as f:
        cPickle.dump(index_end, f)
    with open(save_path + 'position_index.pickle', 'w') as f:
        cPickle.dump(position_index, f)
    print 'data saved'

def create_position_indexes(index_start, index_end):
    """
    create ordered set of gene first NA position
    according to an index per chromosome
    """
    position_index = defaultdict(defaultdict)
    for key in index_start:
        position_index['start'][key] = sorted(index_start[key].keys())
    for key in index_end:
        position_index['end'][key] = sorted(index_end[key].keys())

    return position_index

def load_indexed_gene_annotations(gtf_path=GTF_PATH):
    """
    load index of genes according to chromosomes annotations:
    chr1    unknown exon    11874   12227   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16932";
    chr1    unknown exon    12613   12721   .       +       .       gene_id "DDX11L1"; gene_name "DDX11L1"; transcript_id "NR_046018"; tss_id "TSS16932";
    """
    regex = re.compile('gene\_id "(?P<geneid>.+)"\; gene')
    f = open(gtf_path, "r")
    index_start = defaultdict(dict)
    index_end = defaultdict(dict)

    for line in f:
        line = line.split('\t')

        if int(line[3]) not in index_start[line[0]]:
            index_start[line[0]][int(line[3])] = []

        if int(line[4]) not in index_end[line[0]]:
            index_end[line[0]][int(line[4])] = []


        index_start[line[0]][int(line[3])].append(
            (int(line[4]),
            regex.findall(line[8])[0]))
        index_end[line[0]][int(line[4])].append(
            (int(line[3]),
             regex.findall(line[8])[0]))

    return index_start, index_end


if __name__ =="__main__":
    main()
