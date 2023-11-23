"""
For
- QC on samples.
- split dataset.

Assume fid==iid
"""

# import os
import argparse
from logging import getLogger
from ...system.logger import logger_setting
# import pathlib
import numpy as np
# import pandas as pd
from sklearn.model_selection import ShuffleSplit, KFold
# from sklearn.model_selection import ShuffleSplit, StratifiedShuffleSplit, KFold, StratifiedKFold

from ..io_genot import plink
from ..sample import sampleio as sampleio

# from .ukb import field as ukbfield

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--cv',
                        action='store_true',
                        help='')

    parser.add_argument('--fgenot',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2, plink, etc.')

    parser.add_argument('--seed',
                        type=int,
                        default=None,
                        help='seed of random_state')

    parser.add_argument('--prop-test',
                        type=float,
                        default=None,
                        help='')

    parser.add_argument('--cvn',
                        type=int,
                        default=None,
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def split_samples_multiphe(pheno, prop, seed=None):
    sss = ShuffleSplit(
        n_splits=1, test_size=prop, random_state=seed)
    indexs = list(sss.split(np.zeros_like(pheno.index.to_numpy())))[0]
    logger.debug("indexs[0] {}".format(str(indexs[0][:10])))
    logger.debug("indexs[1] {}".format(str(indexs[1][:10])))
    return indexs


def split_samples_multiphe_cv(sample, cvn, seed=None):
    # pheno = phe['phe'].values
    # print("pheno", pheno)

    skf = KFold(n_splits=cvn, shuffle=True, random_state=seed)
    # list of n_cv tuples, each contains (index_tr),(index_va)
    indexss = list(skf.split(np.zeros_like(sample.index.to_numpy())))
    logger.debug("indexss[0][0] {}".format(str(indexss[0][0])))
    return indexss


def write_cv(args):
    sample = plink.load_sample(args.fgenot, args.genot_format)
    # sample = sample.set_index('id')
    logger.debug('sample\n{}'.format(sample.head()))

    cvn = args.cvn

    # include 'sex' for female-only phenotype
    assert list(sample.columns) == ['id', 'sex']

    cv = sample.set_index('id')

    _, ts = split_samples_multiphe(
        cv, args.prop_test, args.seed)

    # col: 'cv0', ...
    # val: 0: training, 1: validation, 2: test
    for cvi in range(0, cvn):
        # [ref](https://pandas.pydata.org/docs/user_guide/indexing.html#combining-positional-and-label-based-indexing)
        cv_col = 'cv' + str(cvi)
        cv.loc[:, cv_col] = 'NaN'
        cv.iloc[ts, cv.columns.get_loc(cv_col)] = 'ts'

    cv_obs = cv[cv['cv0'] != 'ts']

    indexss_cv = split_samples_multiphe_cv(
        cv_obs, args.cvn, args.seed + 1)
    for cvi, (tr, va) in enumerate(indexss_cv):
        cv_col = 'cv' + str(cvi)
        tr_index = cv_obs.index[tr]
        va_index = cv_obs.index[va]
        cv.loc[tr_index, cv_col] = 'tr'
        cv.loc[va_index, cv_col] = 'va'

    for cvi in range(0, cvn):
        assert len(cv[cv['cv' + str(cvi)] == 'NaN']) == 0

    cv = cv.reset_index(drop=False)
    logger.debug('cv\n{}'.format(cv.head()))
    # cv.to_csv(args.fout, sep='\t', index=False)
    sampleio.write_sample(args.fout, cv)


def main():
    args = argument()

    if args.cv:
        write_cv(args)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
