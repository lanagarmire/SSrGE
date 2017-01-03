from garmire_SSrGE.extract_data import ExtractData
from garmire_SSrGE.load_data import load_gsm_and_sample_names_from_soft

from garmire_SSrGE.config import EXPRESSION_MATRIX_FOLDER_PATH
from garmire_SSrGE.config import VCF_FOLDER_PATH

from sklearn.feature_extraction import DictVectorizer

from os import listdir

from sys import stdout

from tabulate import tabulate

from collections import Counter


def debug():
    """ DEBUG """
    extract_matrix = ExtractMatrix()

    SNV_mat = extract_matrix.extract_SNV_mat()
    GE_mat = extract_matrix.extract_GE_mat()


class ExtractMatrix():
    """
    class to extract SNV_mat and GE_mat from existing dataset

    Project variables must be defined into the config file (config.py):

    PROJECT_PATH
    # path toward the project folder
    GTF_PATH
    # gtf file of the reference genome
    GTF_PATH
    # gtf file of the reference genome
    EXPRESSION_MATRIX_FOLDER_PATH
    # Path of the folders containing the gene expression matrices
    GENE_MATRIX_NAME
    # name of the gene expression matrix file
    VCF_FOLDER_PATH
    # Path of the folders containing the vcf files
    VCF_NAME
    # Name of the VCF the vcf files
    """

    def __init__(self,
                 min_shared_snv=None,
                 min_gene_expr=None,
                 min_average_gene_expr=2,
                 limit=None):
        """
        :min_shared_snv: int    min number of cells sharing a given snv
        :min_gene_expr: float    min number of gene expression value
        :min_average_gene_expr: float    min number of average gene expression value
                                         on average
        """
        self.min_shared_snv = min_shared_snv
        self.min_gene_expr = min_gene_expr
        self.min_average_gene_expr = min_average_gene_expr

        samples_with_vcf = set(listdir(VCF_FOLDER_PATH))
        samples_with_ge_mat = set(listdir(EXPRESSION_MATRIX_FOLDER_PATH))

        self.samples = list(samples_with_vcf.intersection(samples_with_ge_mat))

        if limit:
            self.samples = self.samples[:limit]

        self.samples_snv_dict = {}
        self.samples_ge_dict = {}

        self.extract_data = ExtractData()
        self.gsm_to_name = load_gsm_and_sample_names_from_soft()
        self.names = []

        for sample in self.samples:
            gsm = sample.split('_')[0]
            name = self.gsm_to_name[gsm] if gsm in self.gsm_to_name else gsm
            self.names.append(name)

    def get_samples_list(self):
        """ """
        return self.samples

    def extract_SNV_mat(self):
        """
        construct SNV binary matrix (n_samples x n_SNVs),
        using the project variables described into the config.py file

        return:
            :SNV_mat: Matrix (n_samples x n_SNVs)
        """

        i = 0

        for sample in self.samples:
            self.samples_snv_dict[sample] = self.extract_data.\
                                            load_snv_from_cell(sample)
            i += 1
            stdout.write('\r{0} / {1} VCF files readed'.format(i, len(self.samples)))
            stdout.flush()

        average_snvs = Counter()

        for sample in self.samples_snv_dict:
            average_snvs += self.samples_snv_dict[sample]

        if self.min_shared_snv:
            for sample in self.samples_snv_dict:
                for snv in self.samples_snv_dict[sample].keys():
                    if average_snvs[snv] < self.min_shared_snv:
                        self.samples_snv_dict[sample].pop(snv)

        tab = []

        for sample, name in zip(self.samples, self.names):
            tab.append((sample,
                        name,
                        len(self.samples_snv_dict[sample])))

        print '\n', tabulate(tab, headers=['sample', 'name', 'Number of SNVs'])

        vectorizer = DictVectorizer()

        f_matrix = vectorizer.fit_transform([self.samples_snv_dict[sample]
                                             for sample in self.samples])
        self.snv_index = vectorizer.vocabulary_

        print 'number of SNVs in the dataset:', len(self.snv_index)

        return f_matrix


    def extract_GE_mat(self):
        """
        construct GE matrix (n_genes x n_samples),
        using the project variables described into the config.py file

        return:
            :GE_mat: Matrix  (n_genes x n_samples)
        """

        i = 0

        for sample, name in zip(self.samples, self.names):
            self.samples_ge_dict[sample] = self.extract_data.\
                                           load_expression_profile_from_cell(sample)
            i += 1
            stdout.write('\r{0} / {1} expression files readed'.format(i, len(self.samples)))
            stdout.flush()

        average_expr = self.extract_data.get_average_expression_dict()

        tab = []

        if self.min_gene_expr or self.min_average_gene_expr:
            for sample in self.samples_ge_dict:
                for gene in self.samples_ge_dict[sample].keys():

                    if self.min_average_gene_expr \
                       and average_expr[gene] < self.min_average_gene_expr:
                        self.samples_ge_dict[sample].pop(gene)

                    if self.min_gene_expr \
                       and self.samples_ge_dict[sample][gene] < self.min_gene_expr:
                        self.samples_ge_dict[sample].pop(gene, None)

        for sample, name in zip(self.samples, self.names):
            tab.append((sample,
                        name,
                        len(self.samples_ge_dict[sample])))

        print '\n', tabulate(tab, headers=['sample', 'name', 'Number of genes'])

        vectorizer = DictVectorizer()
        f_matrix = vectorizer.fit_transform([self.samples_ge_dict[sample]
                                             for sample in self.samples])
        self.ge_index = vectorizer.vocabulary_

        print 'number of genes in the dataset:', len(self.ge_index)

        return f_matrix.T


if __name__ == "__main__":
    debug()
