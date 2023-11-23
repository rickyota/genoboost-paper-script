"""
For assoc
"""


import argparse
from logging import getLogger
from ...system.logger import logger_setting
# import decimal
# from decimal import Decimal
import pandas as pd
import numpy as np

from ..io_genot import plink
from ..snv import snvio as snvio
from ..snv import snv as snvop


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        help='output file')

    parser.add_argument('--fassoc',
                        help='')

    parser.add_argument('--fgenot',
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2vzs, plink2, plink, etc.')

    parser.add_argument('--dominance',
                        action='store_true',
                        help='')

    parser.add_argument('--align-allele',
                        action='store_true',
                        help='stay allele of ss')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def add_a2(fout, fassoc):
    def get_a2(x):
        if x['A1'] == x['ALT']:
            return x['REF']
        elif x['A1'] == x['REF']:
            return x['ALT']
        else:
            raise RuntimeError('A1 if not either ALT or REF')

    assoc = snvio.load_assoc(fassoc)

    # remove #
    assoc.columns = assoc.columns.str.replace('#', '')

    if (assoc['ALT'] != assoc['A1']).any():
        raise RuntimeError('Use "omit-ref" option in plink, or comment out this assertion')

    logger.debug('assoc\n{}'.format(assoc.head()))
    assoc.loc[:, 'A2'] = assoc[['REF', 'ALT', 'A1']].apply(lambda x: get_a2(x), axis='columns')

    assoc.to_csv(fout, sep='\t', index=False)


def extract_ss(fout, fassoc, fgenot, genot_format, dominance, align_allele):
    snv = plink.load_snv(fgenot, genot_format)

    assoc = pd.read_csv(fassoc, sep='\t')

    assoc = assoc.rename(columns={'minor_allele': 'a1'})

    vars = assoc['variant'].str.split(':', expand=True)
    assoc.loc[:, 'chrom'] = vars[0]
    assoc.loc[:, 'pos'] = vars[1]
    assoc.loc[:, 'ref'] = vars[2]
    assoc.loc[:, 'alt'] = vars[3]

    def get_a2(x):
        if x['a1'] == x['alt']:
            return x['ref']
        elif x['a1'] == x['ref']:
            return x['alt']
        else:
            raise RuntimeError('A1 if not either ALT or REF', x['a1'], x['alt'], x['ref'], x['variant'])

    assoc.loc[:, 'a2'] = assoc[['variant', 'ref', 'alt', 'a1']].apply(lambda x: get_a2(x), axis='columns')

    assoc = assoc.drop(columns=['variant', 'ref', 'alt'])

    assoc_ext, snv_ext = snvop.match_allele(assoc, snv, both=True, allow_reverse=False, order='left', ignore_index=True,
                                            col1_ref='a2', col1_alt='a1')

    assert (assoc_ext['chrom'] == snv_ext['chrom']).all()
    assert (assoc_ext['pos'] == snv_ext['pos']).all()

    logger.debug('assoc_ext: {}'.format(assoc_ext))

    if dominance:
        col_tstat = 'dominance_tstat'
    else:
        col_tstat = 'tstat'

    # extract tstat only from assoc since it needs reverse
    assoc_ext = assoc_ext.loc[:, ['chrom', 'pos', 'a1', 'a2', col_tstat, 'n_complete_samples']]
    assoc_ext.loc[:, 'id'] = snv_ext.loc[:, 'id']

    if align_allele:
        allele_flip = (assoc_ext['a1'] == snv_ext['ref'])
        logger.debug('flip bet. assoc and plink: {}'.format(np.count_nonzero(allele_flip)))
        assoc_ext.loc[allele_flip, col_tstat] = -assoc_ext.loc[allele_flip, col_tstat]
        assoc_ext.loc[:, 'a1'] = snv_ext.loc[:, 'alt']
        assoc_ext.loc[:, 'a2'] = snv_ext.loc[:, 'ref']

    # rename col for d-ldsc
    assoc_out = assoc_ext.loc[:, ['id', 'a1', 'a2', 'n_complete_samples', col_tstat]]
    assoc_out = assoc_out.rename(columns={'id': 'SNP', 'a1': 'A1', 'a2': 'A2', 'n_complete_samples': 'N', col_tstat: 'Z'})
    assoc_out.to_csv(fout, index=False, sep='\t')


def main():
    args = argument()

    extract_ss(args.fout, args.fassoc, args.fgenot, args.genot_format, args.dominance, args.align_allele)

    # if args.add_a2:
    #    # logger.info('join_snv_hm3]')
    #    add_a2(args.fout, args.fassoc)

    # if args.add_col_p:
    #    # logger.info('join_snv_hm3]')
    #    add_col_p(args.fout, args.fassoc)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
