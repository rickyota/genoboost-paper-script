"""
align 1000g snvs to ukb
"""


import argparse
from logging import getLogger
from ...system.logger import logger_setting
import pandas as pd

from ..io_genot import plink
from ..io_genot import snv as snvio
from ..snv import snv as snvop


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
                        default='plink2vzs',
                        help='plink2, plink, etc.')

    parser.add_argument('--fgenot2',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot2-format',
                        default='plink2vzs',
                        help='plink2vzs, plink2, plink')

    parser.add_argument('--align-snv-rename-1kg-ukb',
                        action='store_true',
                        help='output snv rename file')

    parser.add_argument('--align-snv-rename-sbayesr-ukb',
                        action='store_true',
                        help='output snv rename file')

    parser.add_argument('--fldm',
                        default=None,
                        help='sbayesr')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def align_snv_rename_1kg_ukb(fout, fgenot, genot_format, fgenot2, genot2_format):

    # 1000g
    snv1 = plink.load_snv(fgenot, genot_format)

    # var_rename=snv1[['id']].copy()
    # var_rename.loc[:,'id_new']=''

    var_renames = []

    for chrom in range(1, 23):
        logger.debug('chrom ' + str(chrom))
        # ukb
        snv2 = plink.load_snv(fgenot2, genot2_format, chrom=chrom)

        snv1_extract, snv2_extract = snvop.match_allele(snv1, snv2, order='left', ignore_index=True, both=True, verbose=True)
        logger.debug('snv1\n{}'.format(snv1_extract.head()))
        logger.debug('snv2\n{}'.format(snv2_extract.head()))

        snv_intersect = snv1_extract[['id']].copy().join(snv2_extract[['id']], rsuffix='_new', how='inner')
        assert len(snv_intersect) == len(snv1_extract)
        var_renames.append(snv_intersect)

    var_renames = pd.concat(var_renames, ignore_index=True)

    # 1000g only has rs and ukb has {rs}_{REF}_{ALT} so check if rs are the same
    for rs_old, rs_new in zip(var_renames['id'], var_renames['id_new']):
        if rs_old not in rs_new:
            logger.info('rs_old ' + rs_old)
            logger.info('rs_new ' + rs_new)
            raise RuntimeError("old rs does not contain new rs.")

    var_renames.to_csv(fout, sep='\t', index=False, header=False)


def main():
    args = argument()

    if args.align_snv_rename_1kg_ukb:
        # logger.info('join_snv_hm3]')
        align_snv_rename_1kg_ukb(args.fout, args.fgenot, args.genot_format, args.fgenot2, args.genot2_format)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
