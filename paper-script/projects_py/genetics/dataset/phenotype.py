"""
For phenotype and cov.

Assume fid==iid
"""

# import os
import argparse
from logging import getLogger
from ...system.logger import logger_setting
import pathlib
import numpy as np
import pandas as pd

from ..io_genot import plink
from ..sample import sampleio as sampleio

from .ukb import field as ukbfield

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--cov',
                        action='store_true',
                        help='')

    parser.add_argument('--cov_column',
                        # nargs='+',
                        default=None,
                        help='ex. age,sex,PC_1-PC_10. When use "-", the cov ("PC") should be the same.')

    parser.add_argument('--fgenot',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2, plink, etc.')

    parser.add_argument('--fukb_whole',
                        default=None,
                        help='')

    parser.add_argument('--now_date',
                        nargs='+',
                        help='indicate (year, date)')

    parser.add_argument('--phes',
                        nargs='+',
                        default=None,
                        help='ex. age,sex,PC1-PC10. When use "-", the cov ("PC") should be the same.')

    # phenotype = cov + trait
    # parser.add_argument('--trait',
    #                    action='store_true',
    #                    help='')

    parser.add_argument('--case',
                        action='store_true',
                        help='')

    parser.add_argument('--dukb_case',
                        default=None,
                        help='')

    parser.add_argument('--cov_standard',
                        action='store_true',
                        help='')

    parser.add_argument('--fcov',
                        default=None,
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def convert_ukb_cov(sample, ukb_cov, cov_columns, now_date=None):

    # cov = pd.DataFrame(index=ukb_cov.index)
    cov = sample.set_index('id')

    ukb_cov = ukb_cov.set_index('id')

    cov_names, cov_info = parse_cov_columns_names(cov_columns)

    fieldd = ukbfield.get_cov_fields(cov_names)

    ukb_cov = ukb_cov.loc[cov.index, :].copy()

    for cov_name in cov_names:
        field_cov = fieldd[cov_name]
        info_cov = cov_info[cov_name]

        if cov_name == 'sex':
            assert len(field_cov) == 1
            col = field_cov[0] + '-0.0'
            cov.loc[:, cov_name] = ukb_cov[col]
        elif cov_name == 'age':
            year_now, month_now = tuple([int(x) for x in now_date])
            cov.loc[:, cov_name] = int(year_now) - ukb_cov['34-0.0'].astype(np.int32)
            filt = (ukb_cov['52-0.0'].astype(np.int32) < month_now)
            cov.loc[filt, cov_name] = cov.loc[filt, cov_name] + 1
        elif cov_name == 'PC':
            assert len(field_cov) == 1
            col = field_cov[0]
            start, end = info_cov[0], info_cov[1]
            cov = cov.join(ukb_cov.loc[:, col + '-0.' + start:col + '-0.' + end])

            def rename_col(x, col):
                if x.startswith(col):
                    return 'PC_' + x.split('.')[1]
                else:
                    return x

            # rename col '22009-0.1' -> 'PC_1'
            cov = cov.rename(columns=lambda x: rename_col(x, col))
            # cov = cov.join(ukb_cov.loc[:, 'PC_' + start:'PC_' + end])
    cov = cov.reset_index(drop=False)
    return cov


# TODO: rename 'info'
def parse_cov_columns_names(cov_columns):
    # cov_columns: 'age,PC_1-PC_10'
    # cov_names: ['age','PC']
    # cov_info: {'age': None, 'PC': ('1','10')}
    cov_names = []
    cov_info = {}
    for cov_col in cov_columns.split(','):
        if '-' in cov_col:
            xs = cov_col.split('-')
            assert len(xs) == 2
            x0, x1 = xs[0], xs[1]
            assert '_' in x0
            assert '_' in x1

            cov_name, info0 = x0.split('_')
            cov_name1, info1 = x1.split('_')
            assert cov_name1 == cov_name

            info = (info0, info1)
        else:
            cov_name = cov_col
            info = None
        cov_names += [cov_name]
        cov_info[cov_name] = info
    return cov_names, cov_info


def write_cov(args):
    sample = plink.load_sample(args.fgenot, args.genot_format)
    # sample = sample.set_index('id')
    logger.debug('sample\n{}'.format(sample.head()))

    cov_columns = args.cov_column
    cov_names, _ = parse_cov_columns_names(cov_columns)
    logger.debug('cov_names\n{}'.format(str(cov_names)))

    ukb_cov = ukbfield.load_ukb_cov(args.fukb_whole, cov_names)
    logger.debug('ukb_cov\n{}'.format(ukb_cov.head()))

    cov = convert_ukb_cov(sample, ukb_cov, cov_columns, now_date=args.now_date)
    logger.debug('cov\n{}'.format(cov.head()))

    # cov.to_csv(args.fout, sep='\t', index=False)
    sampleio.write_sample(args.fout, cov)


def load_case(file):
    case = pd.read_csv(file,
                       delim_whitespace=True,
                       header=0,
                       dtype=str,
                       names=['id'])
    return case


def write_case_pheno(args):
    sample = plink.load_sample(args.fgenot, args.genot_format)
    sample = sample[['id']].copy()
    # sample = sample.set_index('id')
    logger.debug('sample\n{}'.format(sample.head()))

    # phe=sample.set_index('id')
    pheno_index = pd.Index(sample['id'])

    pheno = sample.copy()
    pheno = pd.DataFrame(index=pheno_index)
    phes = args.phes

    dcase = pathlib.Path(args.dukb_case)
    for phe in phes:
        # fcase = args.dcase + phe_name + '/case.txt'
        fcase = dcase / phe / 'case.txt'
        case = load_case(fcase)
        # print('case',case['id'])
        common_case = case.set_index('id').index.intersection(pheno_index)
        pheno.loc[:, phe] = '0'
        pheno.loc[common_case, phe] = '1'
        logger.debug("case {}:  {} / {}".format(phe, str(len(pheno[pheno[phe] == '1'])), str(len(pheno))))

    pheno = pheno.reset_index(drop=False)
    # phe.to_csv(args.fout, sep='\t', index=False)
    sampleio.write_sample(args.fout, pheno)


def main():
    args = argument()

    if args.cov:
        # logger.info('qc_sample')
        write_cov(args)

    if args.cov_standard:
        pass

    if args.case:
        write_case_pheno(args)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
