import unittest

from garmire_SSrGE.extract_matrices_from_dataset import ExtractMatrix

from garmire_SSrGE.config import EXPRESSION_MATRIX_FOLDER_PATH
from garmire_SSrGE.config import VCF_FOLDER_PATH
from garmire_SSrGE.config import GTF_PATH
from garmire_SSrGE.config import INDEX_SAVE_PATH

from garmire_SSrGE.generate_refgenome_index import main as make_index

from os.path import isdir
from os.path import isfile

import warnings


class TestPackage(unittest.TestCase):
    """ """
    def test_gtf_path(self):
        """
        test if the GTF path defined in  config exists
        """
        self.assertTrue(isfile(GTF_PATH))

    def test_gtf_index(self):
        """
        test the gtf index creation
        """
        self.assertTrue(make_index())

    def test_gtf_index_path(self):
        """
        test the gtf index path
        """
        self.assertTrue(INDEX_SAVE_PATH)

    def test_vcf_dir_exits(self):
        """
        test if the vcf directory defined in config exists
        """
        self.assertTrue(isdir(VCF_FOLDER_PATH))

    def test_snv_matrices(self):
        """
        test if the snv extraction matrix on one sample
        """
        extract_matrix = ExtractMatrix(limit=1)
        matrix = extract_matrix.extract_SNV_mat()

        if isinstance(matrix, type(None)):
            warnings.warn(
                'SNV matrix is None beacause vcf folder is not defined!')
            return

        self.assertTrue(matrix.shape)

    def test_expression_matrix_dir_exits(self):
        """
        test if the expression matrix directory defined in config exists
        """
        self.assertTrue(isdir(EXPRESSION_MATRIX_FOLDER_PATH))

    def test_ge_matrices(self):
        """
        test if the gene expression extraction matrix on one sample
        """
        extract_matrix = ExtractMatrix(limit=1)
        matrix = extract_matrix.extract_GE_mat()

        if isinstance(matrix, type(None)):
            warnings.warn(
                'gene expression matrix is none beacause GE folder is not defined!')
            return

        self.assertTrue(matrix.shape)


if __name__ == "__main__":
    unittest.main()
