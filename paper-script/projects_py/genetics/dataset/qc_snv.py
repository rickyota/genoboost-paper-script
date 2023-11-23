"""
For
- QC on snps.
"""

import argparse
from logging import getLogger
from ...system.logger import logger_setting
import pandas as pd

from ..io_genot import snv as snvio
from ..snv import snv as snvop


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--fout2',
                        default=None,
                        help='output file')

    parser.add_argument('--fgenot',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2, plink, etc.')

    parser.add_argument('--join_hm3',
                        action='store_true',
                        help='')

    parser.add_argument('--fsnv',
                        default=None,
                        help='')

    parser.add_argument('--fhm3',
                        default=None,
                        help='')

    parser.add_argument('--extract_snvs',
                        action='store_true',
                        help='')

    # indicate fout2
    # parser.add_argument('--join_snvs_both',
    #                    action='store_true',
    #                    help='also output joined snvs on fnsv_join to fout2')

    parser.add_argument('--fsnv_join',
                        default=None,
                        help='assume to have header if has more than one column')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def join_hm3(args):
    # same format as .pvar but no header
    # snv=plink.load_snv_plink2(args.fsnv,header=False)

    # freq = plink.load_freq_plink2(args.ffreq,header=False)
    # logger.debug('freq {}'.format(freq.head() ))

    # use liftover
    hm3 = snvio.load_hm3(args.fhm3)
    logger.debug('hm3 {}'.format(hm3.head()))

    snv = snvio.load_snvs(args.fsnv, header=False)

    snv_intersect = pd.merge(snv, hm3, on=['chrom', 'pos'])

    cols = 'chrom,pos,id,ref,alt'.split(',')
    snv_intersect = snv_intersect[cols]
    snv_intersect.to_csv(args.fout, sep='\t', index=False, header=False)


def extract_snvs_flip_rev_allele(fout1, fsnv1, fsnv2, fout2=None):
    # allow flipped allele
    # fsnv can have header or not with
    # fsnv should at least have ['chrom', 'pos', 'id', 'ref', 'alt']

    snv1 = snvio.load_snvs(fsnv1)
    snv2 = snvio.load_snvs(fsnv2)
    snv1_extract, snv2_extract = snvop.match_allele(snv1, snv2, both=True, verbose=True)
    snv1_extract.to_csv(fout1, sep='\t', index=False, header=False)
    if fout2 is not None:
        snv2_extract.to_csv(fout2, sep='\t', index=False, header=False)


def main():
    args = argument()

    if args.join_hm3:
        # logger.info('join_snv_hm3]')
        join_hm3(args)

    # TODO: which columns to output?
    # if args.join_snvs:
    #    join_snvs_flip_rev_allele(args.fout, args.fsnv, args.fsnv_join)

    if args.extract_snvs:
        extract_snvs_flip_rev_allele(args.fout, args.fsnv, args.fsnv_join, args.fout2)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
