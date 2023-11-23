
import pathlib
import glob
import pandas as pd
import argparse
from logging import getLogger
from ..system.logger import logger_setting

from ..genetics.io_genot import plink
# from ..io import filename as iof
from ..genetics.score import scoreio as scoreio


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--dout',
                        help='output dir')

    parser.add_argument('--dscore',
                        help='')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--phe',
                        help='phenotype name')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


# tmp
def get_fscores(d):
    g_fscore = d / ('ldpred_rho-*.score.gz')
    fscores = glob.glob(str(g_fscore))
    fscores = [pathlib.Path(f) for f in fscores]
    return fscores


def flip_score(dout, dscore, fpheno, phe):
    fscores = get_fscores(dscore)

    for fscore in fscores:

        # score = pd.read_csv(fscore, sep='\t')
        score = scoreio.load_score_fpheno(fscore, fpheno, phe)
        logger.debug('score {}'.format(score.head()))

        # flip for ldpred
        if (score.loc[score['status'] == 1, 'score'].mean()
                < score.loc[score['status'] == 0, 'score'].mean()):

            logger.debug('ldpred score flipped')
            score.loc[:, 'score'] = -score.loc[:, 'score']

        fname = fscore.name
        fout = dout / fname
        score = score.loc[:, ['id', 'score']]
        #print('fout', fout)
        #print('score', score)
        scoreio.write_score(fout, score)


def main():
    args = argument()

    flip_score(pathlib.Path(args.dout), pathlib.Path(args.dscore),
               pathlib.Path(args.fpheno), args.phe)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
