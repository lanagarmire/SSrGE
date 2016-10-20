from garmire_SSrGE.multiprocess_fitting import BatchFitting

from garmire_SSrGE.config import TIME_LIMIT
from garmire_SSrGE.config import MIN_OBS_FOR_REGRESS
from garmire_SSrGE.config import NB_THREADS

from sklearn.linear_model import Lasso

from collections import Counter

from tabulate import tabulate

from scipy.sparse import issparse

import numpy as np


def debug():
    """
    #### DEBUG ####

    **** Test function ****

    """
    from garmire_SSrGE.examples import create_example_matrix_v1

    X, Y, W = create_example_matrix_v1()

    ssrge = SSrGE(alpha=0.01)
    X_r = ssrge.fit_transform(X, Y)
    import ipdb;ipdb.set_trace()
    score = ssrge.score(X,Y)

    ssrge.rank_vSNVs()


class SSrGE():
    """
    Class to perform the SSrGE (Sparse SNV inference to reflect Gene Expression)
    """
    def __init__(
            self,
            time_limit=TIME_LIMIT,
            min_obs_for_regress=MIN_OBS_FOR_REGRESS,
            nb_threads=NB_THREADS,
            model='LASSO',
            model_params=None,
            alpha=0.1,
            l1_ratio=0.5,
            verbose=True):
        """ """
        self.time_limit = time_limit
        self.min_obs_for_regress = min_obs_for_regress
        self.nb_threads = nb_threads

        self.vSNV_weight = None # total vSNV absolute weight
        self.SNV_mat_shape = None # dim of the fitted SNV_mat
        self.GE_mat_shape = None # dim of the GE_mat used as predicat
        self.vSNV_index = None # vSNV index
        self.intercepts = None # Intercepts for non null model: {index gene: intercept value}
        self.coefs = None # coefs for non null model: {index gene: coefs dict}
        self.verbose = verbose # whether to print ranking results into the terminal

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

    def fit(self, SNV_mat, GE_mat, to_dense=False):
        """
        infer vSNV by fitting sparse linear models using SNV as features
        and gene expression as objectives

        input:
            :SNV_mat: (n_samples x n_SNVs) matrix (binary). Matrix can be sparse
            :GE_mat: (n_GE x n_samples) matrix (float value)
            :to_dense: Bool    if True SNV_mat is converted as ndarray

        return:
            SNV_index, vSNV_mat

            :SNV_index: list<int>    List of vSNV index from the SNV_matrix
            :vSNV_mat: (n_samples x n_vSNVs) matrix (binary)    (len(n_vSNVs) < len(n_SNVs))
        """

        self.SNV_mat_shape = SNV_mat.shape
        self.GE_mat_shape = GE_mat.shape

        assert(self.SNV_mat_shape[0] == self.GE_mat_shape[1])

        if issparse(GE_mat):
            GE_mat = GE_mat.todense()

        if to_dense and issparse(SNV_mat):
            SNV_mat = SNV_mat.todense()

        if isinstance(GE_mat, np.matrix):
            GE_mat = np.array(GE_mat)

        g_index, coefs, intercepts = BatchFitting(
            I_mat=SNV_mat,
            O_mat=GE_mat,
            model=self.model,
            model_params=self.model_params,
            nb_processes=self.nb_threads,
            time_limit=self.time_limit,
            min_obs_for_regress=self.min_obs_for_regress,
            only_nonzero=True).run()

        self._process_computed_coefs(coefs, g_index, intercepts)

    def transform(self, SNV_mat):
        """
        create sparse matrix using input original coefs list of dicts
        input:
            :SNV_mat: Matrix (len(samples), len(SNV))
        return:
            :vSNV_mat: Matrix (len(samples), len(vSNV))
        """
        return SNV_mat.T[self.vSNV_weight.keys()].T

    def fit_transform(self, SNV_mat, GE_mat, to_dense=False):
        """
        Combination of fit and transform functions
        """
        self.fit(SNV_mat, GE_mat, to_dense)
        return self.transform(SNV_mat)

    def _process_computed_coefs(self, coefs, g_index, intercepts):
        """
        instanciate weight coefs and vSNV indexes

        input:
            :coefs: list<Counter>
        """
        self.vSNV_index = list(set([key for coef in coefs for key in coef.iterkeys()]))
        self.vSNV_index = {self.vSNV_index[i]: i for i in range(len(self.vSNV_index))}

        self.vSNV_weight = Counter()
        self.intercepts = {}
        self.coefs = {}

        for counter, gene, intercept in zip(coefs, g_index, intercepts):
            self.vSNV_weight += counter

            if counter:
                self.intercepts[gene] = intercept
                self.coefs[gene] = counter

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


            norm = ((Y_inferred - Y_test)**2).sum()
            norm_null = ((Y_null_inferred - Y_test)**2).sum()

            y_n = float(len(Y_inferred))

            errs_model.append(1 / y_n * norm)
            errs_null_model.append(1 / y_n * norm_null)

        return np.mean(errs_model), np.mean(errs_null_model)

    def rank_vSNVs(self, extract_matrix=None):
        """
        rank vSNVs according to their inferred coefs

        input:
            [OPTIONAL] extract_matrix garmire_ssrge.extract_matrix.ExtractMatrix instance
            if passed, allows to parse SNV name and ids to the ranking results
        """

        ranked_snv = sorted(self.vSNV_weight.iteritems(),
                            key=lambda x:x[1],
                            reverse=True)

        snv_ids = ['' for _ in range(len(ranked_snv))]

        names, scores = map(list, zip(*ranked_snv))

        if extract_matrix:
            pos_index = {pos: name for name, pos in extract_matrix.snv_index.iteritems()}

            for i in range(len(names)):
                names[i] = pos_index[names[i]]
                snv_ids[i] = extract_matrix.extract_data.snv_id_dict[names[i]]

        tab = zip(names, snv_ids, scores)

        if self.verbose:
            print '\n'
            print tabulate(tab, headers=['feature name', 'vSNV id', 'score'])

        return tab

    def rank_genes(self, extract_matrix):
        """
        rank genes according to the inferred coefs of vSNVs inferred and present inside

        input:
            extract_matrix garmire_ssrge.extract_matrix.ExtractMatrix instance
        """

        snv_ranked = self.rank_vSNVs(extract_matrix)

        gene_scores = Counter()

        for (gene, pos), snv_id, score in snv_ranked:
            gene_scores[gene] += score

        ranked_genes = sorted(gene_scores.iteritems(),
                              key=lambda x:x[1],
                              reverse=True)

        if self.verbose:
            print tabulate(ranked_genes, headers=['Gene name', 'score'])

        return ranked_genes


if __name__ == "__main__":
    debug()
