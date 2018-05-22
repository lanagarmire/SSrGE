""" """

from sklearn.model_selection import KFold
from garmire_SSrGE.ssrge import SSrGE

from garmire_SSrGE.config import CROSS_VAL_NFOLD

import numpy as np


def debug():
    """
    #### DEBUG ####
    **** Test function ****
    """
    from garmire_SSrGE.examples import create_example_matrix_v1

    X, Y, W = create_example_matrix_v1()

    cross_val = LinearCrossVal(
        model='LASSO',
        SNV_mat=X,
        GE_mat=Y
    )

    path = cross_val.regularization_path('alpha',  [0.01, 0.1, 0.2])

    return path


class LinearCrossVal():
    """
    Class to perform cross-validation
    """
    def __init__(self,
                 SNV_mat,
                 GE_mat,
                 n_folds=CROSS_VAL_NFOLD,
                 verbose=True,
                 **ssrge_params):
        """ """
        if GE_mat.shape[0] == SNV_mat.shape[0] and \
           GE_mat.shape[1] != SNV_mat.shape[1]:
            GE_mat = GE_mat.T

        self.SNV_mat = SNV_mat
        self.GE_mat = GE_mat
        self.verbose = verbose
        self.n_folds = n_folds

        self.ssrge_params = ssrge_params

        self.errs_models = None
        self.errs_empty_models = None
        self.nb_coefs_list = None
        self.nb_models = None
        self.nb_model_mean = None
        self.regularization_value_list = None

    def regularization_path(self, param_name, value_list):
        """
        :param_name: str    the name of the param to test
        :value_list: list(float)
        """
        self.err_model_mean = []
        self.err_empty_mean = []
        self.nb_coefs_mean = []
        self.nb_model_mean = []
        self.intercept_mean = []
        self.regularization_value_list = value_list

        for value in value_list:
            self.ssrge_params[param_name] = value

            (errs_models,
             errs_null_models,
             nb_coefs,
             nb_models,
             intercepts
            ) = self.fit()

            self.err_model_mean.append(errs_models)
            self.err_empty_mean.append(errs_null_models)
            self.nb_model_mean.append(nb_models)
            self.nb_coefs_mean.append(nb_coefs)
            self.intercept_mean.append(intercepts)

            if self.verbose:
                print('\nmean error model:', errs_models)
                print('mean error null model:', errs_null_models)
                print('mean number of model:', nb_models)
                print('mean number of eeSNVs:', nb_coefs)

        return self.err_model_mean

    def fit(self):
        """ """
        i = 0

        errs_models = []
        errs_null_models = []
        nb_coefs = []
        nb_models = []
        intercepts = []

        print('\n######## cross validation\n####parameters:{0}'\
              .format(self.ssrge_params))

        ssrge = SSrGE(**self.ssrge_params)

        kfold = KFold(n_splits=self.n_folds,
                                 shuffle=True)

        for train, test in kfold.split(self.SNV_mat):
            i += 1
            print('\n## fold nb {0}'.format(i))

            X_train = self.SNV_mat[train]
            Y_train = self.GE_mat.T[train].T

            X_test = self.SNV_mat[test]
            Y_test = self.GE_mat.T[test].T

            ssrge.fit(X_train, Y_train)

            score, score_null = ssrge.score(X_test, Y_test)

            errs_models.append(score)
            errs_null_models.append(score_null)
            nb_coefs.append(len(ssrge.eeSNV_weight))
            nb_models.append(len(ssrge.intercepts))
            intercepts.append(np.mean(list(ssrge.intercepts.values())))

        return  (np.mean(errs_models),
                 np.mean(errs_null_models),
                 np.mean(nb_coefs),
                 np.mean(nb_models),
                 np.mean(intercepts)
                 )


if __name__ == '__main__':
    debug()
