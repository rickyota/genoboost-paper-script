
import argparse
import pathlib
import numpy as np
from sklearn import metrics
import pandas as pd
import statsmodels.api as sm
from collections import defaultdict
import glob
from logging import getLogger

from ...system.logger import logger_setting
from ...genetics.sample import sampleio as sampleio
from ...genetics.sample import sample as sampleop
from ..io import filename as iof
from ..score import score as scoreop
from ..wgt import wgt as wgtop

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    # ./result/nonaddpgs/dataset.12/ .. /score.cv0/
    parser.add_argument('--dout',
                        help='output dir')

    parser.add_argument('--method',
                        help='')


    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--phe',
                        help='phe')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def load_fscores(d, method):
    g_fscore = d / (method + '*.score*')
    fscores = glob.glob(str(g_fscore))
    fscores = [pathlib.Path(f) for f in fscores]
    return fscores


def get_score_cols(score):
    # or .startswith('score')
    return [x for x in score.columns if (x != 'id') and (x != 'status')]


def flip_score(dout, fpheno, phe):

    method='ldpred'

    fscores = load_fscores(dout, method)
    logger.debug('fscores {}'.format(fscores))

    for fscore in fscores:
        logger.debug('fscore {}'.format(fscore))
        parad = scoreop.parse_para(fscore)
        logger.debug('parad {}'.format(parad))
        score = scoreop.load_score_phe(fscore, fpheno, phe)
        # print('score', score)
        cols = get_score_cols(score)
        #print('cols', cols)
        # unnecessary to create score_new
        # col='score' or 'score_n100'
        score_new = score[['id']].copy()
        for col in cols:

            # flip for ldpred
            if (score.loc[score['status'] == 1, col].mean()
               < score.loc[score['status'] == 0, col].mean()):
                logger.debug("flipped wgt in col {}".format(col))
                score_new.loc[:, col] = -score[col]
            else:
                score_new.loc[:, col] = score[col]

        score_new.to_csv(fscore.parent / (str(fscore.name) + '_flip'), sep='\t', index=False)


def main():
    args = argument()

    if args.method != 'ldpred':
        raise RuntimeError

    flip_score(pathlib.Path(args.dout),
               args.fpheno,args.phe)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
