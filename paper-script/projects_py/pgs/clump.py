"""
Clump wgt
"""


import argparse
from logging import getLogger
from ..system.logger import logger_setting
import pandas as pd
import numpy as np

from ..genetics.io_genot import plink


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--convert-wgt',
                        action='store_true',
                        help='')

    parser.add_argument('--fclump',
                        default=None,
                        help='')

    parser.add_argument('--fss',
                        default=None,
                        help='')

    parser.add_argument('--ffreq',
                        default=None,
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def load_tsv_header(fin):
    df = pd.read_csv(fin,
                     delim_whitespace=True,
                     header=0,
                     dtype=str)
    return df


def load_clump(fclump):
    clump = load_tsv_header(fclump)
    clump = clump[['CHR', 'SNP', 'BP', 'P']]
    col_map = {'CHR': 'chrom', 'BP': 'pos', 'SNP': 'id', 'P': 'p'}
    clump = clump.rename(columns=col_map)
    return clump


def load_ss(fss):
    ss = load_tsv_header(fss)
    ss = ss[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'A2', 'OR', 'LOG(OR)_SE', 'P']]
    col_map = {'CHROM': 'chrom', 'POS': 'pos', 'ID': 'id', 'REF': 'ref', 'ALT': 'alt',
               'A1': 'a1', 'A2': 'a2', 'OR': 'or', 'LOG(OR)_SE': 'se', 'P': 'p'}
    ss = ss.rename(columns=col_map)
    return ss


def calc_a1_freq(wgt):
    def get_a1_frq(x):
        if x['a1'] == x['alt']:
            return x['alt_frq']
        elif x['a1'] == x['ref']:
            # TODO: format digits
            return str(1.0 - float(x['alt_frq']))
        else:
            raise RuntimeError('A1 if not either ALT or REF')
    wgt.loc[:, 'a1_frq'] = wgt[['a1', 'alt', 'ref', 'alt_frq']].apply(lambda x: get_a1_frq(x), axis='columns')
    return wgt


def convert_wgt(fout, fclump, fss, ffreq):

    clump = load_clump(fclump)
    logger.debug('clump\n{}'.format(clump.head()))
    logger.debug('clump\n{}'.format(clump.info()))

    ss = load_ss(fss)
    logger.debug('ss\n{}'.format(ss.head()))
    logger.debug('ss\n{}'.format(ss.info()))

    ss = ss[['id', 'a1', 'a2', 'or']]

    wgt = pd.merge(clump, ss, on='id', how='left')

    freq = plink.load_freq_plink2(ffreq, header=False)
    freq = freq[['id', 'ref', 'alt', 'alt_frq']]

    wgt = pd.merge(wgt, freq, on='id', how='left')

    wgt = calc_a1_freq(wgt)
    logger.debug('wgt\n{}'.format(wgt.head()))

    wgt.loc[:, 'wgt'] = np.log(wgt.loc[:, 'or'].astype(np.float64))
    wgt = wgt.drop(columns={'or'})

    #print('wgt', wgt.head())

    cols_out = ['chrom', 'pos', 'id', 'a1', 'a2', 'p', 'a1_frq', 'wgt']
    wgt = wgt[cols_out]

    wgt.to_csv(fout, sep='\t', index=False)


def main():
    args = argument()

    if args.convert_wgt:
        logger.info('convert_wgt')
        convert_wgt(args.fout, args.fclump, args.fss, args.ffreq)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
