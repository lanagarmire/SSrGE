import unittest

from os import popen
from os import listdir

from commands import getstatusoutput

from os.path import isfile
from os.path import isdir
from os.path import split as pathsplit

from garmire_SNV_calling.config import OUTPUT_ROOT
from garmire_SNV_calling.config import SOFT_PATH

from garmire_SNV_calling.config import REF_GENOME
from garmire_SNV_calling.config import ANNOTATION_PATH
from garmire_SNV_calling.config import STAR_INDEX_PATH
from garmire_SNV_calling.config import DBSNP
from garmire_SNV_calling.config import VCF_RESOURCES
from garmire_SNV_calling.config import FASTQ_PATH
from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN

from garmire_SNV_calling.config import JAVA
from garmire_SNV_calling.config import GATK_DIR
from garmire_SNV_calling.config import PICARD_DIR
from garmire_SNV_calling.config import PATH_STAR_SOFTWARE

from fnmatch import fnmatch


class TestPackage(unittest.TestCase):
    """ """
    def test_output_root(self):
        """assert that OUTPUT_ROOT folder exists"""
        self.assertTrue(isdir(OUTPUT_ROOT))

    def test_soft(self):
        """assert that .soft file exists"""
        self.assertTrue(isfile(SOFT_PATH))

    def test_ref_genome(self):
        """assert that ref genome file exists"""
        self.assertTrue(isfile(REF_GENOME))

    def test_annotation_path(self):
        """assert that STAR ref folder exists"""
        self.assertTrue(isdir(pathsplit(STAR_INDEX_PATH)[0]))

    def test_dbsnp(self):
        """assert that dbsnp vcf file exists"""
        self.assertTrue(isfile(DBSNP))

    def test_vcf_resources(self):
        """assert that additional vcf file exist"""
        for vcf in VCF_RESOURCES:
            self.assertTrue(isfile(vcf))

    def test_fastq_path(self):
        """assert that fastq path exists"""
        self.assertTrue(isdir(FASTQ_PATH))

    def test_fastq_path_not_empty(self):
        """assert that fastq path exists"""
        self.assertTrue(len(listdir(FASTQ_PATH)))

    def test_fastq_path_with_folders_with_fastqfile(self):
        """assert that fastq folder exists and that .fastq files are inside"""

        i = 0
        for fastq_folder in listdir(FASTQ_PATH):
            if isfile(FASTQ_PATH + fastq_folder):
                continue
            if PATTERN and not fnmatch(fastq_folder, PATTERN):
                continue

            folder = "{0}/{1}".format(FASTQ_PATH, fastq_folder)

            print 'testing if {0} is empty'.format(folder)
            self.assertTrue(filter(lambda fil: fnmatch(fil, '*.fastq'),
                            listdir(folder)))

    def test_java(self):
        """assert that java is installed and > 1.8"""
        res = getstatusoutput('{0} -version'.format(JAVA))[1]
        self.assertIsNotNone(res)

        version = res.split('"')[1].rsplit('.', 1)[0]
        self.assertTrue(float(version) >= 1.8)

    def test_GATK(self):
        """assert that GATK .jar file exists"""
        self.assertTrue(isfile('{0}/GenomeAnalysisTK.jar'.format(GATK_DIR)))

    def test_picard_tools(self):
        """assert that picard-tools .jar files exist"""
        self.assertTrue(isfile('{0}/picard.jar'.format(PICARD_DIR)))
        self.assertTrue(isfile('{0}/picard-lib.jar'.format(PICARD_DIR)))

    def test_STAR_aligner(self):
        """assert that STAR aligner bin exists and return version"""
        self.assertIsNotNone(popen('{0} --version'.format(PATH_STAR_SOFTWARE)))


if __name__ == "__main__":
    unittest.main()
