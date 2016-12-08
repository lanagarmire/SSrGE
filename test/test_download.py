import unittest

from os import popen

from os.path import isfile
from os.path import isdir

from garmire_download_ncbi_sra.config import FASTQ_DUMP
from garmire_download_ncbi_sra.config import PATH_DATA
from garmire_download_ncbi_sra.config import PATH_SOFT

from garmire_download_ncbi_sra.download_data import get_urls

import urllib2

class TestPackage(unittest.TestCase):
    """ """
    def test_fastq_dump(self):
        """assert that fastq-dump exists"""
        self.assertIsNotNone(popen(FASTQ_DUMP))

    def test_is_path(self):
        """assert that data folder exits"""
        self.assertTrue(isdir(PATH_DATA))

    def test_is_soft(self):
        """assert that soft file exits"""
        self.assertTrue(isfile(PATH_SOFT))

    def test_is_urls(self):
        """assert that urls can be extracted from soft files"""
        urls = get_urls()

        self.assertTrue(len(urls))

    def test_connect_to_urls(self):
        """assert that the first url can be reached"""
        urls = get_urls()
        gsm, url = urls[0]

        self.assertTrue(gsm.count('GSM'))
        self.assertTrue(urllib2.urlopen(url))


if __name__ == "__main__":
    unittest.main()
