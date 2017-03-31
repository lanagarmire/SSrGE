#! /usr/bin/python

from os import listdir

from os.path import isdir
from os.path import isfile

from collections import defaultdict
from collections import Counter

from fnmatch import fnmatch
from bisect import bisect
from time import time

import numpy as np

from garmire_SSrGE.config import EXPRESSION_MATRIX_FOLDER_PATH
from garmire_SSrGE.config import GENE_MATRIX_NAME
from garmire_SSrGE.config import VCF_FOLDER_PATH
from garmire_SSrGE.config import VCF_NAME

from garmire_SSrGE.load_data import process_line_from_vcf_file
from garmire_SSrGE.load_data import process_line_from_annotated_vcf_file
from garmire_SSrGE.load_data import load_indexes


def debug():
    """ DEBUG """
    extract_data = ExtractData()


class ExtractData():
    """ """
    def __init__(
            self,
            expression_matrix_folder_path=EXPRESSION_MATRIX_FOLDER_PATH,
            gene_matrix_name=GENE_MATRIX_NAME,
            vcf_folder_path=VCF_FOLDER_PATH,
            vcf_name=VCF_NAME):
        """ """
        self.expression_matrix_folder_path = expression_matrix_folder_path
        self.gene_matrix_name = gene_matrix_name
        self.vcf_folder_path = vcf_folder_path
        self.vcf_name = vcf_name

        self.index = None
        self.position_index = None
        self.index_start = None
        self.index_end = None
        self.snv_id_dict = defaultdict(str)
        self._snvs_index = {}

        self.average_expression = defaultdict(list)

    def _load_indexes(self):
        """ """
        if isinstance(self.index_end, type(None)):
            (self.index_start,
             self.index_end,
             self.position_index) = load_indexes()

    def _load_annotate_snv_from_vcf(self, snv_path):
        """load annotated snv """
        f_snv = open(snv_path, 'r')

        result = defaultdict(list)

        for line in f_snv:
            res = process_line_from_annotated_vcf_file(line)
            if res:
                snvid, snvinfolist = res
                result[snvid] = snvinfolist
        return result

    def load_snv_from_vcf_file(self, vcf_path):
        """
        load snv from a vcf file

        input:
            :vcf_path: path to the vcf file
        """
        self._load_indexes()

        line_nb = sum([1 for l in open(vcf_path, 'r')
                       if l[0] !='#'])

        f_snv = open(vcf_path, 'r')

        wrong_count = 0
        good_count = 0

        result = Counter()

        for line in f_snv:
            res = process_line_from_vcf_file(line)

            if not res:
                continue

            chrid, start, end, snv_id = res

            if chrid not in self.position_index['start']:
                continue

            ref_start = bisect(self.position_index['start'][chrid], start)
            ref_start_list = self.position_index['start'][chrid][max(ref_start-10, 0):
                                                                 ref_start]
            ref_end = bisect(self.position_index['end'][chrid], end)
            ref_end_list = self.position_index['end'][chrid][ref_end:
                                                             ref_end+10]
            ref_start_from_end = set([en[0] for e in ref_end_list
                                      for en in self.index_end[chrid][e]])
            genes_hit_by_snv = ref_start_from_end\
                               .intersection(ref_start_list)

            if genes_hit_by_snv:

                for gene_begin in genes_hit_by_snv:
                    for gene_end_tuple in self.index_start[chrid][gene_begin]:
                        gene_end = gene_end_tuple[0]

                        if not ((gene_begin < start) and (end < gene_end) ):
                            wrong_count += 1
                            continue

                        good_count += 1
                        gene_id = gene_end_tuple[1]
                        snv_name = (gene_id, start)

                        if snv_id:
                            self.snv_id_dict[(gene_id, start)] = snv_id

                        result[snv_name] = 1.0

                        snv_index = (chrid, start)

                        self._snvs_index[snv_index] = snv_name

        return result

    def load_snv_from_cell(self, folder):
        """
        Return SNV found as a dict:
        Counter(snv_id: 1)

        input:
            :folder: str    id of the sample
        """
        f_path = "{0}/{1}/{2}".format(
            self.vcf_folder_path, folder, self.vcf_name)

        return self.load_snv_from_vcf_file(f_path)

    def load_expression_profile_from_cell(self, folder):
        """
        Return cell log FPKM as a dict:
        Counter(gene_id: expr_profile)

        input:
            :folder: str    id of the sample
        """
        f_path = "{0}/{1}/{2}".format(
            self.expression_matrix_folder_path, folder, self.gene_matrix_name)

        fpkm_dict = self.load_expression_profile_from_file(f_path)

        for gid in fpkm_dict:
            fpkm_dict[gid] = np.log(1.0 + fpkm_dict[gid])

        return fpkm_dict

    def get_average_expression_dict(self):
        """ """
        for gid in self.average_expression:

            self.average_expression[gid] = np.mean(
                self.average_expression[gid])

        return self.average_expression

    def load_expression_profile_from_file(self, f_path):
        """
        load expression profile from file
        input:
            :f_path: path to the matrix file
        """
        count_array = defaultdict(list)

        res = Counter()

        f_expr = open(f_path, 'r')
        f_expr.readline()
        f_expr.readline()

        tot_nb_read = 0

        for line in f_expr:
            line = line.strip('\n').split('\t')
            value = float(line[-1])
            g_start = float(line[2].split(';', 1)[0])
            g_end = float(line[3].split(';', 1)[0])
            gene_id = line[0]

            if value:
                tot_nb_read += value
                pre_fpkm = value / (g_end - g_start)
                res[gene_id] = pre_fpkm

            else:
                   count_array[gene_id].append(value)

        for gid in res:
            res[gid] *= 10**9 / tot_nb_read
            count_array[gid].append(res[gid])

            self.average_expression[gid].append(res[gid])

        return res

    def get_snv_id_dict(self):
        """ """
        return defaultdict(str, self.snv_id_dict)


if __name__ == "__main__":
    debug()
