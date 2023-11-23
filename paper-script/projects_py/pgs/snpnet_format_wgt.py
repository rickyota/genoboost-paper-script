

# TODO: create all format_wgt.py for all methods for standart format (icluding order or cols)


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

    parser.add_argument('--dout',
                        help='output dir')
    parser.add_argument('--dwgt',
                        help='')
    parser.add_argument('--fgenot',
                        help='')
    parser.add_argument('--genot-format',
                        default='plink2vzs',
                        help='plink2, plink, etc.')

    parser.add_argument('--sex',
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
def format_wgt(dout, dwgt, fgenot, genot_format, sex):
    # args=sys.argv
    # dout=args[1]

    bim = plink.load_snv(str(fgenot), genot_format)
    bim = bim[['chrom', 'pos', 'id', 'ref', 'alt']].copy()

    # col_map = {'CHROM': 'chrom', 'POS': 'pos', 'ID': 'id', 'REF': 'ref', 'ALT': 'alt'}

    # to load A1:A2
    # if len(sys.argv)==3:
    #    fgt=args[2]
    #    cols='chrom,sida,N1,pos,N2,N3'.split(',')
    #    bim=pd.read_csv(fgt+'.bim',sep='\s+',header=None,names=cols)
    #    print("bim",bim)
    #    #bim.loc[:,'sid']=bim['chrom'].apply(str)+':'+bim['pos'].apply(str)
    #    bim=bim.drop(columns=['chrom','pos','N1','N2','N3'])
    #    print("bim",bim)

    # ls snpnet.n*
    g_fscore = str(dwgt) + '/snpnet.n*.raw'
    fscores = glob.glob(g_fscore)
    print("g_fscore", g_fscore)

    print("fscore", fscores)

    cols = 'var,wgt'.split(',')
    for f in fscores:
        print("fname", f)
        wgt = pd.read_csv(f, sep='\t', header=None, names=cols, dtype=str)

        # print("wgt",wgt)

        if sex == 'female':
            ncov = 11
        elif sex == 'both':
            ncov = 12
        else:
            raise RuntimeError()

        # is_cov=(~wgt['var'].str.contains(':'))

        # wgt_snv=wgt[~is_cov].copy()
        wgt_snv = wgt.iloc[ncov:, :].copy()

        if len(wgt_snv) == 0:
            continue

        # for snvs only
        wgt_snv = create_wgt(wgt_snv, bim)

        # use #snvs not #variant
        nsnvs = len(wgt_snv)

        fout = dout / ('snpnet_n-' + str(nsnvs) + '.wgt')
        # fout=f.replace('.raw','.wgt')
        wgt_snv.to_csv(fout, sep='\t', index=False)

        wgt_cov = wgt.iloc[:ncov, :].copy()
        # wgt_cov=wgt[is_cov].copy()

        # make the format same as .wgtcov
        wgt_cov = wgt_cov.rename(columns={'wgt': 'coef'})
        # wgt_cov=wgt_cov.rename(columns={'wgt':'coef','var':''})

        # wgt_cov=wgt_cov.set_index('var')

        fout_cov = dout / ('snpnet_n-' + str(nsnvs) + '.wgtcov')
        wgt_cov.to_csv(fout_cov, sep='\t', index=False)


def main():
    args = argument()

    logger.info('format_wgt')
    format_wgt(pathlib.Path(args.dout), pathlib.Path(args.dwgt), pathlib.Path(args.fgenot), args.genot_format, args.sex)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
