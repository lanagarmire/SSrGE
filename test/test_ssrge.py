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

        score = ssrge.score(X,Y)

        self.assertTrue(score[0] < score[1])

    def test_ssrge_elasticnet(self):
        """test ssrge procedure with elasticnet model"""
        from garmire_SSrGE.examples import create_example_matrix_v1
        from garmire_SSrGE.ssrge import SSrGE

        X, Y, W = create_example_matrix_v1()
        ssrge = SSrGE(alpha=0.01, model='ElasticNet')

        ssrge.fit(X, Y)
        self.assertTrue(ssrge.eeSNV_weight)

        Xr = ssrge.transform(X)

        self.assertTrue(Xr.sum())
        self.assertTrue(Xr.shape[0] == X.shape[0])
        self.assertTrue(Xr.shape[1] <= X.shape[1])

        snv_ranked = ssrge.rank_eeSNVs()

        self.assertTrue(snv_ranked)

        score = ssrge.score(X,Y)

        self.assertTrue(score[0] < score[1])

    def test_ssrge_cnv(self):
        """test ssrge procedure with cnv matrix"""
        from garmire_SSrGE.examples import create_example_matrix_v3
        from garmire_SSrGE.ssrge import SSrGE

        X, Y, C, W = create_example_matrix_v3()
        ssrge = SSrGE(alpha=0.01)

        ssrge.fit(X, Y, C)
        self.assertTrue(ssrge.eeSNV_weight)

        Xr = ssrge.transform(X)

        self.assertTrue(Xr.sum())
        self.assertTrue(Xr.shape[0] == X.shape[0])
        self.assertTrue(Xr.shape[1] < X.shape[1])

        snv_ranked = ssrge.rank_eeSNVs()

        self.assertTrue(snv_ranked)

        score = ssrge.score(X,Y)

        self.assertTrue(score[0] < score[1])

    def test_ssrge_rank_gene(self):
        """test ssrge and rank genes and snvs"""
        from garmire_SSrGE.examples import create_example_matrix_v2
        from garmire_SSrGE.ssrge import SSrGE

        X, Y, gene_id_list, snv_id_list = create_example_matrix_v2()
        ssrge = SSrGE(
            snv_id_list=snv_id_list,
            gene_id_list=gene_id_list,
            nb_ranked_features=2,
            alpha=0.01)

        ssrge.fit(X, Y)
        self.assertTrue(ssrge.eeSNV_weight)

        Xr = ssrge.transform(X)

        self.assertTrue(Xr.sum())
        self.assertTrue(Xr.shape[0] == X.shape[0])
        self.assertTrue(Xr.shape[1] < X.shape[1])

        snv_ranked = ssrge.rank_eeSNVs()

        self.assertTrue(snv_ranked)

        self.assertTrue(len(ssrge.retained_snvs) == ssrge.nb_ranked_features)
        self.assertTrue(len(ssrge.retained_genes) == ssrge.nb_ranked_features)

        score = ssrge.score(X,Y)

        self.assertTrue(score[0] < score[1])

        subgroup = ssrge.rank_features_for_a_subgroup(range(10))

        self.assertTrue(len(subgroup.gene_expr_distrib[gene_id_list[0]]) == 10)
        self.assertTrue(subgroup.snv_weights_distrib)
        self.assertTrue(subgroup.exp_snv_distrib_comp)

    def test_cross_validation(self):
        """test cross validation procedure"""

        from garmire_SSrGE.linear_cross_validation import debug

        path = debug()
        self.assertTrue(path)


if __name__ == "__main__":
    unittest.main()
