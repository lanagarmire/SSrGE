
import unittest

from os import popen
from os import listdir

from commands import getstatusoutput

from os.path import isfile
from os.path import isdir

from garmire_SNV_calling.config import OUTPUT_ROOT

from garmire_SNV_calling.config import TYPE_VAR
from garmire_SNV_calling.config import FASTQ_PATH
from garmire_SNV_calling.config import SPECIFIC_FILENAME_PATTERN as PATTERN

from garmire_SNV_calling.config import JAVA
from garmire_SNV_calling.config import GATK_DIR
from garmire_SNV_calling.config import GATK_JAR
from garmire_SNV_calling.config import PICARD_DIR
from garmire_SNV_calling.config import PATH_STAR_SOFTWARE
from garmire_SNV_calling.config import SAMTOOLS
from garmire_SNV_calling.config import FREEBAYES

from fnmatch import fnmatch


class TestPackage(unittest.TestCase):
    """ """
    def test_output_root(self):
        """assert that OUTPUT_ROOT folder exists"""
        self.assertTrue(isdir(OUTPUT_ROOT))

    # def test_soft(self):
    #     """assert that .soft file exists"""
    #     self.assertTrue(isfile(SOFT_PATH))

    def test_ref_genome(self):
        """assert that ref genome file exists"""
        for typ in TYPE_VAR:
            self.assertTrue(isfile(TYPE_VAR[typ]['REF_GENOME']))

    # def test_annotation_path(self):
    #     """assert that STAR ref folder exists"""
    #     for typ in TYPE_VAR:
    #         self.assertTrue(isdir(pathsplit(
    #             TYPE_VAR[typ]['STAR_INDEX_PATH'])[0]))

    def test_gtf_file(self):
        """assert that GTF file exists"""
        for typ in TYPE_VAR:
            self.assertTrue(isfile(TYPE_VAR[typ]['ANNOTATION_PATH']))

    def test_dbsnp(self):
        """assert that dbsnp vcf file exists"""
        for typ in TYPE_VAR:
            self.assertTrue(isfile(TYPE_VAR[typ]['DBSNP']))

    def test_vcf_resources(self):
        """assert that additional vcf file exist"""
        for typ in TYPE_VAR:
            for vcf in TYPE_VAR[typ]['VCF_RESOURCES']:
                self.assertTrue(isfile(vcf))

    def test_fastq_path(self):
        """assert that fastq path exists"""
        self.assertTrue(isdir(FASTQ_PATH))

    def test_fastq_path_not_empty(self):
        """assert that fastq path not empty"""
        self.assertTrue(len(listdir(FASTQ_PATH)))

    def test_fastq_path_with_folders_with_fastqfile(self):
        """assert that fastq folder exists and that .fastq files are inside"""

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
        self.assertTrue(isfile('{0}/{1}'.format(GATK_DIR, GATK_JAR)))

    def test_freebayes(self):
        """assert that freebayes file exists"""
        self.assertTrue(isfile(FREEBAYES))

    def test_samtools(self):
        """assert that freebayes file exists"""
        self.assertTrue(isfile(SAMTOOLS))

    def test_picard_tools(self):
        """assert that picard-tools .jar files exist"""
        self.assertTrue(isfile('{0}/picard.jar'.format(PICARD_DIR)))
        self.assertTrue(isfile('{0}/picard-lib.jar'.format(PICARD_DIR)))

    def test_STAR_aligner(self):
        """assert that STAR aligner bin exists and return version"""
        self.assertIsNotNone(popen('{0} --version'.format(PATH_STAR_SOFTWARE)))


if __name__ == "__main__":
    unittest.main()
