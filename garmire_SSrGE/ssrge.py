from garmire_SSrGE.multiprocess_fitting import BatchFitting

from garmire_SSrGE.config import TIME_LIMIT
from garmire_SSrGE.config import MIN_OBS_FOR_REGRESS
from garmire_SSrGE.config import NB_THREADS

from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet

from sklearn.metrics import median_absolute_error

from scipy.stats import fisher_exact

from collections import Counter
from collections import defaultdict

from scipy.sparse import issparse

from warnings import warn

from sys import stdout

import numpy as np


def debug():
    """
    #### DEBUG ####

    **** Test function ****

    """
    from garmire_SSrGE.examples import create_example_matrix_v4

    X, Y, C, ge_list, s_list = create_example_matrix_v4()

    ssrge = SSrGE(snv_id_list=s_list,
                  gene_id_list=ge_list,
                  nb_ranked_features=3,
                  alpha=0.01)
    ssrge.fit_transform(X, Y, C)
    ssrge.score(X, Y)

    print(ssrge.retained_snvs)
    print(ssrge.retained_genes)

    ssrge = SSrGE(nb_ranked_features=2,
                  alpha=0.01)

    ssrge.fit_transform(X, Y, C)

    ssrge.score(X,Y)
    print(ssrge.retained_snvs)
    print(ssrge.retained_genes)


class SSrGE():
    """
    Class to perform the SSrGE (Sparse SNV inference to reflect Gene Expression)
    """
    def __init__(
            self,
            snv_id_list=[],
            gene_id_list=[],
            nb_ranked_features=None,
            time_limit=TIME_LIMIT,
            min_obs_for_regress=MIN_OBS_FOR_REGRESS,
            nb_threads=NB_THREADS,
            model='LASSO',
            model_params=None,
            alpha=0.1,
            l1_ratio=0.5,
            verbose=True,
            **kwargs):
        """
        input:
            :gene_id_list: list of genes ids
            :snv_id_list: list(tuple) <snv ids, gene ids>    the gene ids corresponds
                                                              to the gene where the given
                                                              svn is found
            :nb_ranked_features: int    top ranked features (snvs and genes) to keep
        """
        self.retained_genes = []
        self.retained_snvs = []
        self._do_rank_genes = False
        self._snv_ids_given = False
        self.snv_index = None
        self.gene_index = None
        self.snv_id_dict = None
        self.gene_id_dict = None

        self._cnv_used = None
        self.cnv_score = defaultdict(float)

        self.nb_ranked_features = nb_ranked_features

        if list(snv_id_list):
            try:
                assert(all(len(snv) == 2 for snv in snv_id_list))
            except Exception:
                warn('snv_id_list given but not conform and cannot be used.'\
                     'to rank gene_id_list.'\
                     '\ncorrect format: :snv_id_list: list(tuple) <snv ids, gene ids>')
            else:
                self._do_rank_genes = True
                self._snv_ids_given = True

        self._create_dicts(snv_id_list, gene_id_list)

        self.snvs_ranked = [] # list of tupe (snv, score)
        self.genes_ranked = [] # list of tupe (gene, score)

        self.gene_weights = None

        self.time_limit = time_limit
        self.min_obs_for_regress = min_obs_for_regress
        self.nb_threads = nb_threads

        self.eeSNV_weight = None # total eeSNV absolute weight
        self.SNV_mat = None # fitted SNV_mat
        self.GE_mat = None # fitted GE_mat
        self.SNV_mat_shape = None # dim of the fitted SNV_mat
        self.GE_mat_shape = None # dim of the GE_mat used as predicat
        self.eeSNV_index = None # eeSNV index
        self.intercepts = None # Intercepts for non null model: {index gene: intercept value}
        self.coefs = None # coefs for non null model: {index gene: coefs dict}
        self.verbose = verbose # whether to print ranking results into the terminal
        self.eeSNV_CIS_score = defaultdict(float)
        self.gene_CIS_score = defaultdict(float)

        if model == 'LASSO':
            self.model = Lasso
            self.model_params = {
                'alpha': alpha,
                'max_iter': 1000,
            }
        elif model == 'ElasticNet':
            self.model = ElasticNet
            self.model_params = {
                'alpha': alpha,
                'l1_ratio': l1_ratio,
                'max_iter': 1000,
            }
        else:
            self.model = model
            self.model_params = model_params

    def _create_dicts(self, snv_id_list, gene_id_list):
        """ """
        self.snv_index = dict(enumerate(snv_id_list))
        self.gene_index = dict(enumerate(gene_id_list))

        self.snv_id_dict = {name: pos
                            for pos, name in self.snv_index.iteritems()}
        self.gene_id_dict = defaultdict(str, {name: pos
                            for pos, name in self.gene_index.iteritems()})

    def fit(self, SNV_mat, GE_mat, CNV_mat=None, to_dense=False):
        """
        infer eeSNV by fitting sparse linear models using SNV as features
        and gene expression as objectives

        input:
            :SNV_mat: (n_samples x n_SNVs) matrix (binary). Matrix can be sparse
            :GE_mat: (n_GE x n_samples) matrix (float value)
            :to_dense: Bool    if True SNV_mat is converted as ndarray

        return:
            SNV_index, eeSNV_mat

            :SNV_index: list<int>    List of eeSNV index from the SNV_matrix
            :eeSNV_mat: (n_samples x n_eeSNVs) matrix (binary)    (len(n_eeSNVs) < len(n_SNVs))
        """

        self.SNV_mat_shape = SNV_mat.shape
        self.GE_mat_shape = GE_mat.shape

        if not self._snv_ids_given:
            self._create_dicts(range(self.SNV_mat_shape[1]),
                               range(self.GE_mat_shape[0]))

        assert(self.SNV_mat_shape[0] == self.GE_mat_shape[1])
        assert(self.SNV_mat_shape[0] == self.GE_mat_shape[1])

        if issparse(GE_mat):
            GE_mat = GE_mat.todense()

        if to_dense and issparse(SNV_mat):
            SNV_mat = SNV_mat.todense()

        if isinstance(GE_mat, np.matrix):
            GE_mat = np.array(GE_mat)

        self._cnv_used = CNV_mat is not None

        g_index, coefs, intercepts = BatchFitting(
            I_mat=SNV_mat,
            O_mat=GE_mat,
            CNV_mat=CNV_mat,
            model=self.model,
            model_params=self.model_params,
            nb_processes=self.nb_threads,
            time_limit=self.time_limit,
            min_obs_for_regress=self.min_obs_for_regress,
            only_nonzero=True).run()

        self._process_computed_coefs(coefs, g_index, intercepts)
        self._rank_eeSNVs()

        if self._do_rank_genes:
            self._rank_genes()

        self.select_top_ranked_features()

        self.SNV_mat = SNV_mat
        self.GE_mat = GE_mat

    def select_top_ranked_features(self, nb_ranked_features=None):
        """ """
        if not nb_ranked_features:
            nb_ranked_features=self.nb_ranked_features

        self.retained_genes = [gene for gene, score in
                               self.genes_ranked[:nb_ranked_features]]
        self.retained_snvs = [snv for snv, score in
                               self.snvs_ranked[:nb_ranked_features]]

    def transform(self, SNV_mat):
        """
        create sparse matrix using input original coefs list of dicts
        input:
            :SNV_mat: Matrix (len(samples), len(SNV))
        return:
            :eeSNV_mat: Matrix (len(samples), len(eeSNV))
        """
        return SNV_mat.T[[self.snv_id_dict[snv] for snv in self.retained_snvs]].T

    def fit_transform(self, SNV_mat, GE_mat, CNV_mat=None, to_dense=False):
        """
        Combination of fit and transform functions
        """
        self.fit(SNV_mat, GE_mat, CNV_mat=CNV_mat, to_dense=to_dense)
        return self.transform(SNV_mat)

    def _process_computed_coefs(self, coefs, g_index, intercepts):
        """
        instanciate weight coefs and eeSNV indexes

        input:
            :coefs: list<Counter>
        """
        if self.verbose:
            print('\nprocess computed coefs....')

        self.eeSNV_index = list(set([key for coef in coefs for key in coef.iterkeys()]))
        self.eeSNV_index = {self.eeSNV_index[i]: i for i in range(len(self.eeSNV_index))}

        self.eeSNV_weight = defaultdict(float)
        self.eeSNV_CIS_score = defaultdict(float)

        self.intercepts = {}
        self.coefs = {}

        i = 0
        length = len(coefs)

        for counter, gene, intercept in zip(coefs, g_index, intercepts):
            gene_name = self.gene_index[gene] if self.gene_index else gene

            for key, count in counter.iteritems():
                if self._cnv_used and key == self.SNV_mat_shape[1]:
                    self.cnv_score[gene_name] = count
                    continue

                self.eeSNV_weight[key] += count

                if self._snv_ids_given:
                    genename, pos = self.snv_index[key]

                    if genename == self.gene_index[gene]:
                        self.eeSNV_CIS_score[self.snv_index[key]] += count

            if counter:
                self.intercepts[gene] = intercept
                self.coefs[gene] = counter

                if self.cnv_score[gene_name]:
                    self.coefs[gene].pop(self.SNV_mat_shape[1])

            i += 1

            stdout.write('\r {0:.2f} / 100'.format(i / length * 100))
            stdout.flush()

        for snv in self.eeSNV_CIS_score:
            self.eeSNV_CIS_score[snv] /= self.eeSNV_weight[self.snv_id_dict[snv]]

        print('\n')

    def score(self, SNV_mat, GE_mat):
        """
        Return mean of MSE for GE prediction, using GE_mat as predicat
        and SNV_mat as features. SNV_mat and GE_mat should be havethe same number
        of SNVs and genes than the fitted models, respectively

        input:
            :SNV_mat: (n_samples x n_SNVs) matrix (binary). Matrix can be sparse
            :GE_mat: (n_GE x n_samples) matrix (float value)

        return:
            :err_models: float    mean of the MSE for models
            :err_null_models: float    mean of the MSE for null models (only intercepts)
        """
        assert(SNV_mat.shape[1] == self.SNV_mat_shape[1])
        assert(GE_mat.shape[0] == self.GE_mat_shape[0])

        errs_model = []
        errs_null_model = []

        if issparse(GE_mat):
            GE_mat = GE_mat.todense()

        if isinstance(GE_mat, np.matrix):
            GE_mat = np.array(GE_mat)

        if issparse(SNV_mat):
            SNV_mat = SNV_mat.todense()

        if isinstance(SNV_mat, np.matrix):
            SNV_mat = np.array(SNV_mat)

        for non_null_gene in self.coefs:
            non_zero = np.nonzero(GE_mat[non_null_gene])[0]

            if not len(non_zero):
                continue

            Y_test = GE_mat[non_null_gene][non_zero]
            coef = np.zeros(self.SNV_mat_shape[1])
            coef[self.coefs[non_null_gene].keys()] = [self.coefs[non_null_gene][k]
                                                      for k in self.coefs[non_null_gene]]
            Y_inferred =np.asarray(SNV_mat[non_zero] * np.matrix(coef).T).T[0] \
                        + self.intercepts[non_null_gene]

            Y_null_inferred = np.ones(Y_test.shape[0]) * self.intercepts[non_null_gene]

            score = median_absolute_error(Y_inferred, Y_test)
            score_null = median_absolute_error(Y_null_inferred, Y_test)

            errs_model.append(score)
            errs_null_model.append(score_null)

        return np.mean(errs_model), np.mean(errs_null_model)

    def rank_eeSNVs(self):
        """
        rank eeSNVs according to their inferred coefs
        """
        return self.snvs_ranked

    def _rank_eeSNVs(self):
        """
        rank eeSNVs according to their inferred coefs
        """

        self.snvs_ranked = []

        ranked_snv = sorted(self.eeSNV_weight.iteritems(),
                            key=lambda x:x[1],
                            reverse=True)

        for snv_i, score in ranked_snv:
            self.snvs_ranked.append((self.snv_index[snv_i], score))

        return self.snvs_ranked

    def rank_genes(self):
        """
        rank genes according to their inferred coefs
        """
        return self.genes_ranked

    def _rank_genes(self):
        """
        rank genes according to the inferred coefs of eeSNVs inferred and present inside
        """
        self.gene_weights = defaultdict(float)
        self.gene_CIS_score = defaultdict(float)

        for snv_i, score in self.eeSNV_weight.iteritems():
            snv = self.snv_index[snv_i]
            gene, pos = snv
            self.gene_weights[gene] += score

            if self._snv_ids_given:
                self.gene_CIS_score[gene] += self.eeSNV_CIS_score[snv] * score

        for gene in self.gene_CIS_score:
            self.gene_CIS_score[gene] /= self.gene_weights[gene]

        self.genes_ranked = sorted(self.gene_weights.iteritems(),
                                  key=lambda x:x[1],
                                  reverse=True)
        return self.genes_ranked

    def rank_features_for_a_subgroup(self, sample_id_list):
        """
        Rank the eeSNVs and the genes for a given subgroup of samples

        input:
            :sample_id_list: id of samples of interest
                             example [1,5,10] => group with samples 1, 5 and 10
        output:
            :SubGroupData: data container with features specific to the subgroup
        """
        gene_weights_list = defaultdict(Counter)
        snv_weights_list = defaultdict(Counter)
        exp_gene_weights_list = defaultdict(Counter)
        exp_snv_weights_list = defaultdict(Counter)

        sample_id_comp = list(set(
            range(self.SNV_mat.shape[0])).difference(sample_id_list))

        SNV_mat_sub = self.SNV_mat[sample_id_list].todense()
        SNV_mat_comp = self.SNV_mat[sample_id_comp].todense()
        GE_mat_sub = self.GE_mat.T[sample_id_list]
        GE_mat_comp = self.GE_mat.T[sample_id_comp]

        for key in self.eeSNV_weight:
            SNV_mat_sub.T[key] *= self.eeSNV_weight[key]

        subgroup = SubGroupData()

        for index, gene in self.gene_index.iteritems():
            subgroup.gene_expr_distrib[gene] = GE_mat_sub.T[index]

        for snv_i, score in self.eeSNV_weight.iteritems():
            snv = self.snv_index[snv_i]
            gene, pos = snv
            gene_i = self.gene_id_dict[gene]

            for cell_i in xrange(SNV_mat_sub.shape[0]):
                snv_weights_list[cell_i][snv] = SNV_mat_sub[cell_i, snv_i]

                if gene_i != '':
                    gene_weights_list[cell_i][gene] += SNV_mat_sub[cell_i, snv_i]

                if gene_i != '' and GE_mat_sub[cell_i, gene_i]:
                    exp_snv_weights_list[cell_i][snv] = SNV_mat_sub[cell_i, snv_i]
                    exp_gene_weights_list[cell_i][gene] += SNV_mat_sub[cell_i, snv_i]

            if gene_i:
                index_cells_comp = np.nonzero(GE_mat_comp.T[gene_i])[0]
                subgroup.exp_snv_distrib_comp[snv] = np.array(
                    SNV_mat_comp.T[snv_i, index_cells_comp])

        for cell_i in snv_weights_list:

            for gene in gene_weights_list[cell_i]:
                subgroup.gene_weights_distrib[gene].append(
                    gene_weights_list[cell_i][gene])

            for snv in snv_weights_list[cell_i]:
                subgroup.snv_weights_distrib[snv].append(
                    snv_weights_list[cell_i][snv])

        for cell_i in exp_snv_weights_list:
            for gene in exp_gene_weights_list[cell_i]:
                subgroup.exp_gene_weights_distrib[gene].append(
                    exp_gene_weights_list[cell_i][gene])

            for snv in exp_snv_weights_list[cell_i]:
                subgroup.exp_snv_weights_distrib[snv].append(
                    exp_snv_weights_list[cell_i][snv])

        subgroup._get_significant_subgroup_features()

        return subgroup

class SubGroupData():
    """
    class containing data for a given subgroup of cells

    attribute:
        :significant_eeSNVs: list of ranked significant eeSNVs with their score
        :significant_genes: list of ranked significant eeSNVs with their score

        :gene_expr_distrib: distribution of gene expression for each gene
        :gene_weights_distrib: distribution of gene weights for each gene
                               (according to the eeSNVs) for the subgroup
        :snv_weights_distrib: distribution of the eeSNV weights for each eeSNV
        :exp_gene_weights_distrib: distribution of gene weights for each gene using only,
                                   for a given gene, the subset of cells expressing the gene
        :exp_snv_weights_distrib: distribution of eeSNV weights for each eeSNV using only,
                                   for a given eeSNV, the subset of cells expressing the gene
                                   related to the eeSNV
    """
    def __init__(self):
        """ """
        self.significant_eeSNVs = []
        self.significant_genes = []
        self.ranked_eeSNVs = []
        self.ranked_genes = []

        self.gene_expr_distrib = defaultdict(list)
        self.gene_weights_distrib = defaultdict(list)
        self.snv_weights_distrib = defaultdict(list)
        self.exp_gene_weights_distrib = defaultdict(list)
        self.exp_snv_weights_distrib = defaultdict(list)
        self.exp_snv_distrib_comp = defaultdict(list)


    def _get_significant_subgroup_features(self, thres=0.05):
        """ """
        snv_ranked = []
        gene_ranked = defaultdict(float)

        for snv in self.snv_weights_distrib:
            distrib_test = np.asarray(self.exp_snv_weights_distrib[snv]).astype('bool')
            distrib_ref = np.asarray(self.exp_snv_distrib_comp[snv]).astype('bool')

            key_mean = np.mean(self.exp_snv_weights_distrib[snv])
            contingency = np.array([[distrib_test.sum(), (distrib_test == False).sum()],
                                    [distrib_ref.sum(), (distrib_ref == False).sum()],])

            score, pvalue = fisher_exact(contingency)

            if pvalue < thres and distrib_test.mean() > distrib_ref.mean():
                snv_ranked.append((snv, key_mean, pvalue))

        snv_ranked.sort(key=lambda x:x[1], reverse=True)

        self.significant_eeSNVs = snv_ranked

        for (gene, pos), score, pvalue in snv_ranked:
            gene_ranked[gene] += score

        self.significant_genes = sorted(gene_ranked.items(),
                                        key=lambda x:x[1],
                                        reverse=True)

        self.ranked_genes = sorted([(gene, np.mean(self.gene_weights_distrib[gene]))
                                    for gene in self.gene_weights_distrib],
                                   key=lambda x:x[1],
                                   reverse=True)
        self.ranked_eeSNVs = sorted([(snv, np.mean(self.snv_weights_distrib[snv]))
                                     for snv in self.snv_weights_distrib],
                                    key=lambda x:x[1],
                                    reverse=True)


if __name__ == "__main__":
    debug()
