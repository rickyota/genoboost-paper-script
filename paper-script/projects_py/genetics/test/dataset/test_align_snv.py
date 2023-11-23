

import unittest
from io import StringIO
import pandas as pd
import numpy as np


from ...snv import snv as snvop


class TestDataset(unittest.TestCase):

    def test_match_allele(self):
        snv1 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'A', 'C'),
            ('1', '2', 'A', 'C'),
            ('1', '3', 'A', 'C'),
            ('1', '4', 'A', 'C'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        snv2 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'C', 'A'),
            ('1', '2', 'T', 'G'),
            ('1', '3', 'G', 'T'),
            ('1', '4', 'A', 'T'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        ans1 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'A', 'C'),
            ('1', '2', 'A', 'C'),
            ('1', '3', 'A', 'C'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        ans2 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'C', 'A'),
            ('1', '2', 'T', 'G'),
            ('1', '3', 'G', 'T'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        ans_filt1 = np.array([True, True, True, True, False])
        ans_filt2 = np.array([True, True, True, True, False])

        df1, filt1, df2, filt2 = snvop.match_allele(snv1, snv2, both=True, return_filter=True)

        print('df1', df1)
        print('df2', df2)
        print('ans_filt1', ans_filt1)
        print('ans_filt2', ans_filt2)

        pd.testing.assert_frame_equal(ans1, df1)
        pd.testing.assert_frame_equal(ans2, df2)

        np.testing.assert_array_equal(ans_filt1, filt1)
        np.testing.assert_array_equal(ans_filt2, filt2)

    def test_match_allele_order(self):
        snv1 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'A', 'C'),
            ('1', '2', 'A', 'C'),
            ('1', '3', 'A', 'C'),
            ('1', '4', 'A', 'C'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        snv2 = pd.DataFrame([
            ('1', '3', 'G', 'T'),
            ('1', '1', 'C', 'A'),
            ('1', '2', 'T', 'G'),
            ('1', '4', 'A', 'T'),
            ('1', '0', 'A', 'C'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        ans1 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'A', 'C'),
            ('1', '2', 'A', 'C'),
            ('1', '3', 'A', 'C'),
        ], columns=['chrom', 'pos', 'ref', 'alt'])

        ans2 = pd.DataFrame([
            ('1', '0', 'A', 'C'),
            ('1', '1', 'C', 'A'),
            ('1', '2', 'T', 'G'),
            ('1', '3', 'G', 'T'),
        ], columns=['chrom', 'pos', 'ref', 'alt'], index=[4, 1, 2, 0])

        ans_filt1 = np.array([True, True, True, True, False])
        ans_filt2 = np.array([True, True, True, False, True])

        df1, filt1, df2, filt2 = snvop.match_allele(snv1, snv2, both=True, return_filter=True, order='left')

        print('df1', df1)
        print('df2', df2)
        print('ans_filt1', ans_filt1)
        print('ans_filt2', ans_filt2)

        pd.testing.assert_frame_equal(ans1, df1)
        pd.testing.assert_frame_equal(ans2, df2)

        np.testing.assert_array_equal(ans_filt1, filt1)
        np.testing.assert_array_equal(ans_filt2, filt2)
        # pd.testing.assert_frame_equal(ans_filt1, filt1)
        # pd.testing.assert_frame_equal(ans_filt2, filt2)


if __name__ == '__main__':
    unittest.main()
