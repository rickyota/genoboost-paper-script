"""
For
- conversion of sida to rs
"""

import argparse
from logging import getLogger
from ...system.logger import logger_setting
import pandas as pd

from ..io_genot import plink
# from ..io_genot import snv as snvio
# from ..snv import snv as snvop


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--fgenot',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2, plink, etc.')
    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def sida_rs(fout, fgenot, genot_format):

    snv = plink.load_snv(fgenot, genot_format)

    logger.debug('snv\n{}'.format(snv.head()))

    # sida=chr:pos:ref:alt
    snv1 = snv.copy()
    # sida=chr:pos:alt:ref
    snv2 = snv.copy()

    snv1 = snv1.assign(sida=snv1.apply(lambda x: '%s:%s:%s:%s' % (x['chrom'], x['pos'], x['ref'], x['alt']), axis=1))
    snv2 = snv2.assign(sida=snv2.apply(lambda x: '%s:%s:%s:%s' % (x['chrom'], x['pos'], x['alt'], x['ref']), axis=1))

    df_rs = pd.concat([snv1, snv2])
    # df_rs = pd.concat([snv1, snv2], ignore_index=True)

    df_rs = df_rs.sort_index()
    # df_rs = df_rs.sort_values(['chrom', 'pos'])

    # df_rs = df_rs.drop(columns=['chrom', 'pos', 'ref', 'alt'])
    df_rs = df_rs[['sida', 'id']].copy()

    df_rs = df_rs.reset_index(drop=True)

    logger.debug('df_rs\n{}'.format(df_rs.head()))

    df_rs.to_csv(fout, sep='\t', index=False)


def main():
    args = argument()

    sida_rs(args.fout, args.fgenot, args.genot_format)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
