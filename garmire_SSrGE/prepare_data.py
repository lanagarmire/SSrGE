#! /usr/bin/python

from os import listdir
from os import mkdir
from os import popen

from os.path import isdir
from os.path import isfile
from time import time

from copy import copy
from fnmatch import fnmatch
from collections import Counter
from collections import defaultdict

from sys import stdout

import numpy as np
import re
import cPickle

from scipy.sparse import vstack
from scipy.sparse import csc_matrix

from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction.text import TfidfTransformer

from tabulate import tabulate

from compute_eSNV_distance.extract_data\
    import ExtractData

from compute_eSNV_distance.load_data import load_folder_to_id_dict
from compute_eSNV_distance.load_data import load_aligner_unique_read
from compute_eSNV_distance.load_data import load_deduplicated_check
from compute_eSNV_distance.load_data import load_filter_patterns

from compute_eSNV_distance.config import PATH_DATA_SNV
from compute_eSNV_distance.config import PATH_DATA_ALIGNER
from compute_eSNV_distance.config import PATH_INDEX
from compute_eSNV_distance.config import MIN_SNV_FEATURE
from compute_eSNV_distance.config import MIN_EXPR_FEATURE
from compute_eSNV_distance.config import MIN_GENE_EXPRESSION
from compute_eSNV_distance.config import MIN_AVERAGE_GENE_EXPR
from compute_eSNV_distance.config import MIXED_COEF
from compute_eSNV_distance.config import LIMIT
from compute_eSNV_distance.config import TF_IDF
from compute_eSNV_distance.config import MIX_WITH_EXPR
from compute_eSNV_distance.config import CHECK_DEDUPLICATED
from compute_eSNV_distance.config import CHECK_UNIQUE_READS
from compute_eSNV_distance.config import MIN_SHARED_NB_SNV
from compute_eSNV_distance.config import LOAD_SNV_ANNOTATION
from compute_eSNV_distance.config import PICKLELIZE
from compute_eSNV_distance.config import REMOVE_OLD_PICKLE
from compute_eSNV_distance.config import LOAD_SNV_READ_COVERAGE
from compute_eSNV_distance.config import SAVE_RANKING

from compute_eSNV_distance.config import PATH_DATA_TMP


from compute_eSNV_distance.cluster_labels import ClusterLabels


class PrepareData():
    """class to extract data from snv anf GE files """
    def __init__(self, print_tab=True):
        self.expression_profile_analysis = None
        self.cell_list =  None
        self.gsm_to_id = {}
        self.expr_dicts = {}
        self.cell_dicts = {}
        self.snv_coverage_dicts = {}
        self.nonzero_snvs = Counter()
        self.nonzero_genes = Counter()
        self.annotated_snv_dicts = defaultdict(list)
        self.cosmic_id_dict = defaultdict(str)
        self.average_cell_dict = Counter()
        self.cluster_labels = ClusterLabels()
        self.average_expr_dict = None
        self._load_pickle = PICKLELIZE
        self._rm_old_pickle = REMOVE_OLD_PICKLE

        self.snv_matrix_position = None
        self.snv_index = {}
        self.ge_matrix_position = None
        self.gene_index = {}
        self.snv_id_dict = defaultdict(str)
        self.print_tab = print_tab
        self._tab_printed = False
        self.cell_id_dict = {}
        self.id_list = []
        self.cell_ids = []
        self.f_matrix_snv = None
        self.f_matrix_ge = None
        self._data_removed_already = False

        self._mktmpdir()
        self._prepare_feature_extraction()

    def _prepare_feature_extraction(self):
        """ """
        self.gsm_to_id = load_folder_to_id_dict()
        self.cell_list =  [cell for cell in listdir(PATH_DATA_SNV)
                           if cell.split('__')[0] in self.gsm_to_id]

        self._remove_cell_ids()

        self._perform_sample_QC_check()
        self.extracted_data = ExtractData()

        if self._rm_old_pickle:
            self._rm_old_pickle_dicts()

        if not self._load_pickle_dicts():
            self.extracted_data.load_indexes()
            self.extracted_data.load_cosmic_db()
            self._create_dicts()
            self.extracted_data._cosmic_db = {}

        print 'number of annotated SNVs:', len(self.snv_id_dict)
        print 'number of COSMIC SNVs:', len(self.cosmic_id_dict)

    def _create_dicts(self):
        """ """
        i = 0

        for cell in self.cell_list:
            i += 1

            self.expr_dicts[cell] = self.extracted_data\
                        .load_expression_profile_from_cell(cell)

            f_path_snv = "{0}/{1}/snv_filtered.vcf"\
                         .format(PATH_DATA_SNV, cell)
            f_path_snv_annotated = "{0}/{1}/snv_filtered_annotated.vcf"\
                                   .format(PATH_DATA_SNV, cell)

            self.cell_dicts[cell] = self.extracted_data\
                                        .load_snv_from_vcf(f_path_snv)

            if LOAD_SNV_ANNOTATION:
                self.annotated_snv_dicts[cell] = self.extracted_data\
                                                     .load_annotate_snv_from_vcf(
                                                         f_path_snv_annotated)

            stdout.write('\r {0} / {1} cells loaded'.format(i, len(self.cell_list)))
            stdout.flush()

        vectorizer = DictVectorizer()
        f_matrix = vectorizer.fit_transform(self.cell_dicts.values())

        f_matrix = f_matrix.T.tocsr()

        for snv in vectorizer.vocabulary_:
            self.average_cell_dict[snv] = f_matrix[
                vectorizer.vocabulary_[snv]].sum()
        print '\ndict update done...'

        self.average_expr_dict = self.extracted_data\
                            .get_average_expression_dict()
        print 'expression profile loaded...'

        self.snv_id_dict = self.extracted_data.get_snv_id_dict()
        self.cosmic_id_dict = self.extracted_data.get_cosmic_id_dict()

        self._save_pickle_dicts()

    def _get_snv_read_coverage(self, snvs_list):
        """
        for each cell get the read coverage
        of each SNVs retained in the feature matrix
        """
        if self._load_snv_depth_pickle_dicts():
            return

        print 'Read depth for each SNVs in the dataset...'

        self.extracted_data._make_bed_file(snvs_list)
        nb_cell = len(self.cell_list)
        i = 0
        t = time()

        for cell in self.cell_list:
            f_path_bam = "{0}/{1}/Aligned.sortedByCoord.out.bam"\
                         .format(PATH_DATA_ALIGNER, cell)
            self.snv_coverage_dicts[cell] = self.extracted_data\
                                                .load_snv_read_coverage(
                                                    f_path_bam)
            i += 1
            stdout.write('\r {0} / {1} cells done'.format(i, nb_cell))
            stdout.flush()

        print 'done in {0} s'.format(time() - t)
        self._save_snv_depth_dicts()

    def _perform_sample_QC_check(self):
        """
        perform quality control on fastq files and bam files
        according to aligner output stats previously computed
        two files are checked in PATH_OUTPUT dir:
            - aligner_unique_read.csv
            - deduplicated_check.csv
        """
        to_remove = []

        if CHECK_DEDUPLICATED:
            quality_dict = load_deduplicated_check()

            for folder in self.cell_list:
                if quality_dict[folder] == 'fail':
                    to_remove.append(folder)

        if CHECK_UNIQUE_READS:
            quality_dict = load_aligner_unique_read()

            for folder in self.cell_list:
                if quality_dict[folder] < CHECK_UNIQUE_READS:
                    to_remove.append(folder)

        self.cell_list = set(self.cell_list).difference(to_remove)

        if LIMIT:
            self.cell_list = set(list(self.cell_list)[:LIMIT])


    def _remove_cell_ids(self):
        """
        remove cell ids from cell_list if doesnt match pattern
        """
        to_keep, to_remove = [], []
        filter = load_filter_patterns()

        if filter['TO_KEEP']:
            for pattern in filter['TO_KEEP']:
                for cell in self.cell_list:
                    gsm = cell.split('__', 1)[0]
                    if re.findall(pattern, self.gsm_to_id[gsm]):
                        to_keep.append(cell)

            self.cell_list = set(to_keep)

        for pattern in filter['TO_REMOVE']:
            for cell in self.cell_list:
                gsm = cell.split('__', 1)[0]
                if re.findall(pattern, self.gsm_to_id[gsm]):
                    to_remove.append(cell)

        self.cell_list = set(self.cell_list).difference(to_remove)

    def snv_extraction(self):
        """
        extract features from SVC file output from eSNV
        param:
            :use_high_conf_only: bool   wether to include low coverage SNV additionaly
        return:
            id_list, d_matrix
        """
        feature_dicts = []
        id_list = []
        tab = []

        i = 0

        for cell in self.cell_list:
            i += 1

            cell_gsm = cell.split('__')[0]
            cell_dict = copy(self.cell_dicts[cell])
            expr_dict = self.expr_dicts[cell]

            for snv_id in cell_dict.keys():
                gid, position = snv_id

                if self.nonzero_snvs and not self.nonzero_snvs[snv_id]:
                    cell_dict.pop(snv_id)
                elif gid and (not expr_dict[gid] or \
                (expr_dict[gid] < MIN_GENE_EXPRESSION \
                   or self.average_expr_dict[gid] < MIN_AVERAGE_GENE_EXPR)):
                    cell_dict.pop(snv_id)
                elif MIN_SHARED_NB_SNV and \
                     self.average_cell_dict[snv_id] < MIN_SHARED_NB_SNV:
                    cell_dict.pop(snv_id)
                else:
                    if MIX_WITH_EXPR:
                        cell_dict[snv_id] *= expr_dict[gid]

            if MIN_SNV_FEATURE and len(cell_dict) < MIN_SNV_FEATURE:
                continue

            cell_id = self.gsm_to_id[cell_gsm]
            id_list.append(cell)
            tab.append((cell_id, cell, len(cell_dict)))
            feature_dicts.append(copy(cell_dict))

        if self.print_tab and not self._tab_printed:
            print tabulate(tab, headers=["sample name", "cell id", "feature nb"])
            self._tab_printed = True

        print "\n ====> number of cells:{0}".format(len(tab))

        vectorizer = DictVectorizer()
        f_matrix = vectorizer.fit_transform(feature_dicts)

        self.snv_matrix_position = vectorizer.vocabulary_
        print 'number of SNVs in the dataset:', len(self.snv_matrix_position)
        self.snv_index = {j:i for i,j in self.snv_matrix_position.iteritems()}

        if TF_IDF:
            tfidf = TfidfTransformer()
            f_matrix = tfidf.fit_transform(f_matrix)

        if LOAD_SNV_READ_COVERAGE:
            self._get_snv_read_coverage(self.snv_matrix_position)

        self._mk_cell_id_dict(id_list)
        self.f_matrix_snv = f_matrix

        return id_list, f_matrix

    def gene_expression_extraction(self):
        """
        extract expression features from SVC file output from eSNV
        :use_high_conf_only: bool    wether to include low coverage SNV additionaly
        NO TF-IDF performed!
        return:
            id_list, d_matrix
        """
        feature_dicts = []
        id_list = []
        tab = []
        i = 0

        for cell in self.cell_list:
            i += 1

            cell_dict = copy(self.expr_dicts[cell])
            cell_gsm = cell.split('__')[0]

            for gid in cell_dict.keys():
                if self.nonzero_genes and not self.nonzero_genes[gid]:
                    cell_dict.pop(gid)
                elif cell_dict[gid] < MIN_GENE_EXPRESSION \
                   or self.average_expr_dict[gid] < MIN_AVERAGE_GENE_EXPR:
                    cell_dict.pop(gid)

            if MIN_EXPR_FEATURE and len(cell_dict) < MIN_EXPR_FEATURE:
                continue

            cell_id = self.gsm_to_id[cell_gsm]

            id_list.append(cell)
            tab.append((cell_id, len(cell_dict)))
            feature_dicts.append(copy(cell_dict))

        if self.print_tab and not self._tab_printed:
            print tabulate(tab, headers=["file", "feature nb"])
            self._tab_printed = True

        print "\n ====> number of cells:{0}".format(len(tab))

        vectorizer = DictVectorizer()
        f_matrix = vectorizer.fit_transform(feature_dicts)

        self.ge_matrix_position = vectorizer.vocabulary_
        print 'number of genes in the dataset:', len(self.ge_matrix_position)
        self.gene_index = {j:i for i,j in self.ge_matrix_position.iteritems()}

        if TF_IDF:
            tfidf = TfidfTransformer()
            f_matrix = tfidf.fit_transform(f_matrix)

        self._mk_cell_id_dict(id_list)
        self.f_matrix_ge = f_matrix

        return id_list, f_matrix

    def extract_common_feature_matrix(self):
        """
        Extract feature matrix from both SNV and GE
        for common cells presents in both dataset
        """
        #### extract features ####
        id_list_snv, feature_matrix_snv = self.snv_extraction()
        id_list_ge, feature_matrix_ge = self\
                                          .gene_expression_extraction()

        id_intersection = list(set(id_list_snv).intersection(id_list_ge))

        index_snv = {id: i for i, id in enumerate(id_list_snv)}
        index_ge = {id: i for i, id in enumerate(id_list_ge)}

        #### use only common cell ####
        feature_matrix_r_snv = vstack([feature_matrix_snv[index_snv[id]]
                                       for id in id_intersection])
        feature_matrix_r_ge = vstack([feature_matrix_ge[index_ge[id]]
                                       for id in id_intersection])

        self._mk_cell_id_dict(id_intersection)
        self.f_matrix_snv = feature_matrix_r_snv
        self.f_matrix_ge = feature_matrix_r_ge

        return self.f_matrix_snv, self.f_matrix_ge, id_intersection

    def get_cell_ids(self, id_list):
        """ return experiment id according to GSM"""
        cell_id = [self.gsm_to_id[folder.split('__', 1)[0]]
                                  for folder in id_list]
        return cell_id

    def _mk_cell_id_dict(self, id_list):
        """ """
        self.id_list = id_list
        self.cell_id_dict = {j:i for i,j in enumerate(id_list)}

    def _mktmpdir(self):
        """ """
        if not isdir(PATH_DATA_TMP):
            mkdir(PATH_DATA_TMP)

    def _load_pickle_dicts(self):
        """ """
        if not self._load_pickle or \
           not isfile(PATH_DATA_TMP + 'cell_dicts.pickle'):
            return

        t = time()
        self.cell_dicts = cPickle.load(open(PATH_DATA_TMP + 'cell_dicts.pickle'))
        self.expr_dicts = cPickle.load(open(PATH_DATA_TMP + 'expr_dicts.pickle'))
        self.average_cell_dict = cPickle.load(open(PATH_DATA_TMP +
                                                   'average_cell_dict.pickle'))
        self.annotated_snv_dicts = cPickle.load(open(PATH_DATA_TMP +
                                                     'annotated_snv_dicts.pickle'))
        self.average_expr_dict = cPickle.load(open(PATH_DATA_TMP +
                                                   'average_expr_dict.pickle'))
        self.snv_id_dict = defaultdict(str, cPickle.load(open(PATH_DATA_TMP +
                                                              'snv_id_dict.pickle')))
        self.cosmic_id_dict = defaultdict(str, cPickle.load(open(PATH_DATA_TMP +
                                                   'cosmic_id_dict.pickle')))

        print 'data loaded from pickle files in {0} s...'.format(time() - t)

        return True

    def _save_pickle_dicts(self):
        """ """
        if not self._load_pickle:
            return

        t = time()
        cPickle.dump(self.cell_dicts,
                     open(PATH_DATA_TMP + 'cell_dicts.pickle', 'w'))
        cPickle.dump(self.expr_dicts,
                     open(PATH_DATA_TMP + 'expr_dicts.pickle', 'w'))
        cPickle.dump(self.average_cell_dict,
                     open(PATH_DATA_TMP + 'average_cell_dict.pickle', 'w'))
        cPickle.dump(self.annotated_snv_dicts,
                     open(PATH_DATA_TMP + 'annotated_snv_dicts.pickle', 'w'))
        cPickle.dump(self.average_expr_dict,
                     open(PATH_DATA_TMP + 'average_expr_dict.pickle', 'w'))
        cPickle.dump(self.snv_id_dict,
                     open(PATH_DATA_TMP + 'snv_id_dict.pickle', 'w'))
        cPickle.dump(self.cosmic_id_dict,
                     open(PATH_DATA_TMP + 'cosmic_id_dict.pickle', 'w'))
        print 'data picklized in {0} s...'.format(time() - t)

    def _save_snv_depth_dicts(self):
        """ """
        if not self._load_pickle:
            return

        cPickle.dump(self.snv_coverage_dicts,
                     open(PATH_DATA_TMP + 'snv_coverage_dicts.pickle', 'w'))

    def _load_snv_depth_pickle_dicts(self):
        """ """
        if not self._load_pickle or \
           not isfile(PATH_DATA_TMP + 'snv_coverage_dicts.pickle'):
            return

        self.snv_coverage_dicts = cPickle.load(open(PATH_DATA_TMP +
                                                   'snv_coverage_dicts.pickle'))
        print 'SNVs depth loaded'
        return True

    def _load_sparse_linear_coefs(self, plot_name):
        """ """
        coefs = None

        plot_name = "{0}_coefs.npz".format(plot_name)
        f_path = "{0}/npz/{1}".format(PATH_DATA_TMP, plot_name)

        if isfile(f_path) and SAVE_RANKING:
            print 'loading SNVs coefs from npz file...'
            t = time()
            loader = np.load(f_path)
            coefs = csc_matrix(
                (loader['data'],
                 loader['indices'],
                 loader['indptr']),
                shape=loader['shape'])
            print 'done in {0} s'.format(time() - t)

        return coefs

    def _save_sparse_linear_coefs(self, coefs, plot_name):
        """ """
        if not SAVE_RANKING:
            return

        if not isdir("{0}/npz/".format(PATH_DATA_TMP)):
            mkdir("{0}/npz/".format(PATH_DATA_TMP))

        plot_name = "{0}_coefs".format(plot_name)
        f_path = "{0}/npz/{1}".format(PATH_DATA_TMP, plot_name)

        t = time()
        np.savez(f_path,
                 data=coefs.data,
                 indices=coefs.indices,
                 indptr=coefs.indptr,
                 shape=coefs.shape)
        print '\nsparse linear coefs saved in {0} s'.format(time() - t)

    def _rm_old_ranking(self):
        """ """
        if isdir("{0}/npz/".format(PATH_DATA_TMP)):
            print 'removing npz files...'
            cmd = 'rm {0}/npz/*.npz'.format(PATH_DATA_TMP)
            print popen(cmd).read()

    def _rm_old_pickle_dicts(self):
        """ """
        if self._data_removed_already:
            return

        print 'removing pickle files...'
        cmd = 'rm {0}/*.pickle'.format(PATH_DATA_TMP)
        print popen(cmd).read()
        self._data_removed_already = True

    def _set_nonzero_counters_empty(self):
        """ """
        self.nonzero_snvs = Counter()
        self.nonzero_genes = Counter()

    def _generate_true_labels(self):
        """ """
        self.cell_ids = self.get_cell_ids(self.id_list)
        self.true_labels = self.cluster_labels.generate_cluster_id(self.cell_ids)
