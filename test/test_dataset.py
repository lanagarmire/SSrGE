import unittest


class TestDataset(unittest.TestCase):
    """ """
    def test_dataset(self):
        """ Test dataset pipeline procedure """
        from garmire_SSrGE.examples import launch_pipeline_and_rank_genes

        ranked_genes = launch_pipeline_and_rank_genes(limit=10)

        self.assertTrue(ranked_genes)
