""" example """

from scipy.sparse import csr_matrix
import numpy as np


def create_example_matrix_v1(nb_cells=100, nb_snvs=6, nb_genes=5):
    """
    create a random feature matrix and infer Y according to coefs W
    Four sparse coefs are set into W
    """
    X = csr_matrix(np.random.random((nb_cells, nb_snvs)))
    W = np.zeros((nb_snvs, nb_genes))


    W[0][0] = 5
    W[1][0] = 5
    W[3][3] = 2
    W[5][4] = 6

    Y = (X * W).T

    return X, Y, W

def create_example_matrix_v2(nb_cells=100, nb_snvs=6, nb_genes=5):
    """
    create a random feature matrix and infer Y according to coefs W
    Four sparse coefs are set into W
    create fake snv list and fake gene list
    """
    X = csr_matrix(np.random.random((nb_cells, nb_snvs)))
    W = np.zeros((nb_snvs, nb_genes))

    gene_id_list = range(nb_genes)

    snv_id_list = [(i, i) if i < nb_genes else (0, i)
                   for i in range(nb_snvs)]

    W[0][0] = 5
    W[1][0] = 5
    W[3][3] = 2
    W[5][4] = 6

    Y = (X * W).T

    return X, Y, gene_id_list, snv_id_list

def launch_pipeline_and_rank_genes(alpha=0.1, limit=None):
    """
    launch the SSrGE pipeline on the dataset defined in config.py
    and rank genes and eeSNVs

    input:
        :alpha: float    regularization parameter
    """
    from garmire_SSrGE.extract_matrices_from_dataset import ExtractMatrix
    from garmire_SSrGE.ssrge import SSrGE

    extract_matrix = ExtractMatrix(min_shared_snv=3, limit=limit)
    SNV_mat = extract_matrix.extract_SNV_mat()
    GE_mat = extract_matrix.extract_GE_mat()

    ssrge = SSrGE(alpha=alpha, min_obs_for_regress=5)

    ssrge.fit(SNV_mat, GE_mat)

    ranked_genes = ssrge.rank_genes(extract_matrix)

    score, score_null = ssrge.score(SNV_mat, GE_mat)

    print '\nerror model:', score
    print 'error null model:', score_null

    return ranked_genes


if __name__ == '__main__':
    launch_pipeline_and_rank_genes()
