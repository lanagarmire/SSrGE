import pysam

from os import mkdir
from os.path import isdir
from os.path import isfile

from sys import stdout

from os import popen

from collections import defaultdict

import cPickle

from multiprocessing import Pool

from glob import glob

from os import popen

from os import remove


######################## VARIABLE ############################
PATH_FASTQ = '/mnt/nas_rna2/opoirion/10x_data/fastq/'
PATH_BAM =  '/mnt/nas_rna2/opoirion/10x_data/neurons_900_possorted_genome_bam.bam'

# maximum of reads for a cells
MAX_READS = None
# number of cells you want
NB_CELLS = 1000
FASTQ_THREAD = 10
###############################################################


def main():
    """
    """
    stats, cell_list = get_cell_stats()
    write_bam_files(cell_list)
    bam_to_fastq(cell_list)


def bam_to_fastq(cell_list):
    """
    """
    print('converting all the bam files into fastq files...')
    pool = Pool(FASTQ_THREAD)
    pool.map(_bam_to_fastq, cell_list)

def _bam_to_fastq(cell):
    """
    """
    bam_files = glob('{0}/{1}/*.bam'.format(PATH_FASTQ, cell))

    if not bam_files:
        return

    for bam_file in  bam_files:
        cmd = 'bamToFastq -i {0} -fq {1}/{2}/{2}.fastq'.format(bam_file, PATH_FASTQ, cell)
        popen(cmd).read()
        remove(bam_file)

def write_bam_files(cell_list):
    """
    """
    print('\n#### SECOND PASS ####')
    cmd = "samtools idxstats {0} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}'".format(PATH_BAM)
    nb_reads = int(popen(cmd).read().strip('\n'))

    if MAX_READS is None:
        max_reads = nb_reads
    else:
        max_reads = MAX_READS

    file_dict = {}

    f_raw = pysam.AlignmentFile(PATH_BAM, 'rb')
    header = f_raw.header

    i = 0

    while i < max_reads:
        try:
            reads = f_raw.next()
        except Exception:
            break

        try:
            bc_tags = reads.get_tag('CB')
        except KeyError:
            continue

        i += 1

        stdout.write('\r nb reads {0} / {1}'.format(i, nb_reads))
        stdout.flush()

        if bc_tags not in cell_list:
            continue

        if bc_tags not in file_dict:
            folder = '{0}/{1}'.format(PATH_FASTQ, bc_tags)

            if not isdir(folder):
                mkdir(folder)

            file_dict[bc_tags] = pysam.AlignmentFile(
                '{0}/bc_tags.bam'.format(folder), 'wb', header=header)

        file_dict[bc_tags].write(reads)

    f_raw.close()

def get_cell_stats():
    """
    """
    path_pickle = '{0}/cell_stats.pickle'.format(PATH_FASTQ)

    if isfile(path_pickle):
        stats_dict = cPickle.load(open(path_pickle))
    else:
        stats_dict = _get_cell_stats()
        cPickle.dump(stats_dict, open(path_pickle, 'w'))

    cells, count = zip(*sorted(stats_dict.items(), key=lambda x:x[1], reverse=True)[:NB_CELLS])

    return stats_dict, set(cells)

def _get_cell_stats():
    """
    """
    print('#### FIRST PASS ####')
    cmd = "samtools idxstats {0} | awk -F '\t' '{{s+=$3+$4}}END{{print s}}'".format(PATH_BAM)
    nb_reads = int(popen(cmd).read().strip('\n'))

    if MAX_READS is None:
        max_reads = nb_reads
    else:
        max_reads = MAX_READS

    stats_dict = defaultdict(int)

    f_raw = pysam.AlignmentFile(PATH_BAM, 'rb')

    i = 0

    while i < max_reads:
        try:
            reads = f_raw.next()
        except Exception:
            break

        try:
            bc_tags = reads.get_tag('CB')
        except KeyError:
            continue

        stats_dict[bc_tags] += 1

        i += 1

        stdout.write('\r nb reads {0} / {1}'.format(i, nb_reads))
        stdout.flush()

    f_raw.close()

    return stats_dict


if __name__ == '__main__':
    main()
