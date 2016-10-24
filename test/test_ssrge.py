import unittest


class TestPackage(unittest.TestCase):
    """ """
    def test_multiprocess(self):
        """ Test multiprocess procedure """
        from garmire_SSrGE.multiprocess_fitting import debug

        g_index, coefs, intercepts = debug()

        self.assertTrue(g_index)
        self.assertTrue(coefs)
        self.assertTrue(coefs[0])
        self.assertTrue(sum(intercepts))


    def test_ssrge(self):
        """test ssrge procedure"""
        from garmire_SSrGE.examples import create_example_matrix_v1
        from garmire_SSrGE.ssrge import SSrGE

        X, Y, W = create_example_matrix_v1()
        ssrge = SSrGE(alpha=0.01)

        ssrge.fit(X, Y)
        self.assertTrue(ssrge.eeSNV_weight)

        Xr = ssrge.transform(X)

        self.assertTrue(Xr.sum())
        self.assertTrue(Xr.shape[0] == X.shape[0])
        self.assertTrue(Xr.shape[1] < X.shape[1])

        snv_ranked = ssrge.rank_eeSNVs()

        self.assertTrue(snv_ranked)


    def test_cross_validation(self):
        """test cross validation procedure"""

        from garmire_SSrGE.linear_cross_validation import debug

        path = debug()
        self.assertTrue(path)


if __name__ == "__main__":
    unittest.main()
