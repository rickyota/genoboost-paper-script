
import pathlib
import glob
import pandas as pd
import argparse
from logging import getLogger
from ..system.logger import logger_setting

from ..genetics.io_genot import plink


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        help='output dir')
    parser.add_argument('--fwgt',
                        help='')
    parser.add_argument('--fconvert',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def create_wgt(wgt, bim):

    def func_sida(v):
        vs = v.split('_')
        return pd.Series(vs, index=['id', 'a1_lasso'])

    wgt.loc[:, ['id', 'a1_lasso']] = wgt['var'].apply(func_sida)
    # wgt.loc[:,['sida','A1_lasso']]=wgt['var'].apply(func_sida)

    wgt = pd.merge(wgt, bim, how='left', on='id')

    # check if 'sida' contains A1:A2
    # vid_split_n=len(wgt.iloc[0].loc['vid'].split(':'))
    # if vid_split_n==4:
    #    wgt.loc[:,'sida']=wgt['vid']
    # elif vid_split_n==2:
    #    # add 'sida'
    #    wgt=pd.merge(wgt,bim,how='left',left_on='vid',right_on='sid')
    #    #print("wgt",wgt)
    # else:
    #    raise RuntimeError("sth wrong")

    # create chrom,pos,A1,A2,wgt
    # def func_split(v):
    #    vs=v.split(':')
    #    #print("vs",vs)
    #    return pd.Series(vs,index=['chrom','pos','A1','A2'])

    # wgt.loc[:,['chrom','pos','A1','A2']]=wgt['sida'].apply(func_split)

    # print("wgt",wgt)

    flip = (wgt['a1_lasso'] == wgt['ref'])
    assert (wgt.loc[~flip, 'a1_lasso'] == wgt.loc[~flip, 'alt']).all(), "why?"
    wgt.loc[flip, 'wgt'] = -wgt.loc[flip, 'wgt']
    wgt.loc[:, 'a1'] = wgt['alt']
    wgt.loc[:, 'a2'] = wgt['ref']

    wgt = wgt[['id', 'chrom', 'pos', 'a1', 'a2', 'wgt']].copy()

    # flip sign if allele is filpped
    # flip=(wgt['a1_lasso']==wgt['a2'])
    # assert (wgt.loc[~flip,'a1']==wgt.loc[~flip,'Aa_lasso']).all(), "why?"

    # wgt.loc[flip,'wgt']=-wgt.loc[flip,'wgt']

    # wgt=wgt.drop(columns=['var','a1_lasso'])

    # print("wgt",wgt)

    # align cols
    # wgt=wgt[['sida','chrom','pos','A1','A2','wgt']]

    return wgt


# format wgt to input classic_score
def format_wgt(fout, fwgt, fconvert):

    convert = pd.read_csv(fconvert, sep='\t')

    logger.debug('convert\n{}'.format(convert.head()))

    sida_to_rs = dict(zip(convert['sida'], convert['id']))

    wgt = pd.read_csv(fwgt, sep='\t', dtype=str)

    wgt = wgt.assign(sida=wgt.apply(lambda x: '%s:%s:%s:%s' % (x['chrom'], x['pos'], x['A1'], x['A2']), axis=1))

    wgt.loc[:, 'id'] = wgt['sida'].map(sida_to_rs)

    # remove 'sid' for rust score
    wgt = wgt.drop(columns=['sid', 'sida'])

    logger.debug('wgt\n{}'.format(wgt.head()))

    wgt.to_csv(fout, sep='\t', index=False)


def main():
    args = argument()

    logger.info('format_wgt')
    format_wgt(pathlib.Path(args.fout), pathlib.Path(args.fwgt), pathlib.Path(args.fconvert))


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
