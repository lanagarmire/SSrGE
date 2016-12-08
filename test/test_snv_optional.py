import unittest

from garmire_SNV_calling.config import JAVA
from garmire_SNV_calling.config import SNPEFF
from garmire_SNV_calling.config import FEATURE_COUNT
from garmire_SNV_calling.config import FASTQC

from commands import getstatusoutput


class TestPackage(unittest.TestCase):
    """ """
    def test_snpeff(self):
        """assert that snpEff is installed"""
        self.assertFalse(getstatusoutput("{0} -jar {1} -version".format(JAVA, SNPEFF))[0])

    def test_featurecount(self):
        """assert that featureCounts is installed"""
        self.assertFalse(getstatusoutput("{0} -v".format(FEATURE_COUNT))[0])

    def test_fastqc(self):
        """assert that fastQC is installed"""
        self.assertFalse(getstatusoutput("{0} -version".format(FASTQC))[0])


if __name__ == "__main__":
    unittest.main()
