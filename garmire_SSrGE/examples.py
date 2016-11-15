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


if __name__ == '__main__':
    launch_pipeline_and_rank_genes()
