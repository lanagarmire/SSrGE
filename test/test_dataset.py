import unittest

from garmire_SSrGE.config import GLOBAL_DATA_ROOT
from garmire_SSrGE.config import PROJECT_PATH
from garmire_SSrGE.config import SOFT_PATH
from garmire_SSrGE.config import GTF_PATH
from garmire_SSrGE.config import GENE_MATRIX_NAME
from garmire_SSrGE.config import EXPRESSION_MATRIX_FOLDER_PATH
from garmire_SSrGE.config import VCF_FOLDER_PATH
from garmire_SSrGE.config import VCF_NAME

from os.path import isfile
from os.path import isdir

from os import listdir


class TestDataset(unittest.TestCase):
    """ """
    def test_output_root(self):
        """assert that GLOBAL_DATA_ROOT folder exists"""
        self.assertTrue(isdir(GLOBAL_DATA_ROOT))

    def test_project_path(self):
        """assert that PROJECT_PATH folder exists"""
        self.assertTrue(isdir(PROJECT_PATH))

    def test_gtf_path(self):
        """assert that GTF_PATH file exists"""
        self.assertTrue(isfile(GTF_PATH))

    def test_GE_path(self):
        """assert that EXPRESSION_MATRIX_FOLDER_PATH folder exists"""
        self.assertTrue(isdir(EXPRESSION_MATRIX_FOLDER_PATH))

    def test_GE_files(self):
        """assert that Gene expression files exits"""
        for folder in listdir(EXPRESSION_MATRIX_FOLDER_PATH):
            f_path = EXPRESSION_MATRIX_FOLDER_PATH + folder

            print 'testing that {0} is a folder...'.format(f_path)
            self.assertTrue(isdir(f_path))

            ge_file = '{0}/{1}'.format(f_path, GENE_MATRIX_NAME)
            print 'testing that {0} is a file...'.format(ge_file)
            self.assertTrue(isfile(ge_file))

    def test_VCF_path(self):
        """assert that VCF_FOLDER_PATH folder exists"""
        self.assertTrue(isdir(VCF_FOLDER_PATH))

    def test_VCF_files(self):
        """assert that VCF files exits"""
        for folder in listdir(VCF_FOLDER_PATH):
            f_path = VCF_FOLDER_PATH + folder

            print 'testing that {0} is a folder...'.format(f_path)
            self.assertTrue(isdir(f_path))

            vcf_file = '{0}/{1}'.format(f_path, VCF_NAME)
            print 'testing that {0} is a file...'.format(vcf_file)
            self.assertTrue(isfile(vcf_file))

    def test_soft(self):
        """assert that .soft file exists"""
        self.assertTrue(isfile(SOFT_PATH))


    def test_dataset(self):
        """ Test dataset pipeline procedure """
        from garmire_SSrGE.examples import launch_pipeline_and_rank_genes

        ranked_genes = launch_pipeline_and_rank_genes(limit=10)

        self.assertTrue(ranked_genes)
