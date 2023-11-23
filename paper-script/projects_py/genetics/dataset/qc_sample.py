"""
For
- QC on samples.
- split dataset.

Assume fid==iid
"""

import argparse
from logging import getLogger
from ...system.logger import logger_setting
import numpy as np
import pandas as pd

from ..io_genot import plink

from .ukb import field as ukbfield

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--qc_sample',
                        action='store_true',
                        help='qc on samples')

    parser.add_argument('--fgenot',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2vzs, plink2, plink, etc.')

    parser.add_argument('--fukb_whole',
                        default=None,
                        help='ukb whole dataset')

    parser.add_argument('--fwithdrawn',
                        default=None,
                        help='')

    parser.add_argument('--ethnic',
                        default=None,
                        help='None=eur / noneur')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def write_qc_sample(args):
    # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    # fam = fam.set_index('iid')
    # fam = fam.sort_index()

    if args.ethnic is None:
        ethnic = 'eur'
    else:
        ethnic = args.ethnic

    sample = plink.load_sample(args.fgenot, args.genot_format, raw=True)
    sample = sample.set_index('id')

    # sex is NaN for invalid iid
    # print(sample['sex'].unique())
    # print(sample.loc[~((sample['sex']==1)|(sample['sex']==2)),:])

    logger.debug('sample\n{}'.format(sample.head()))

    logger.info('bfr all qc {}'.format(len(sample)))

    sample = extract_valid_id(sample)
    logger.info('afr valid id {}'.format(len(sample)))

    sample = sample.astype({'sex': np.int32})
    logger.debug('sample\n{}'.format(sample.head()))

    # TODO
    # withdrawn = pd.read_csv(args.fwithdrawn)
    sample = exclude_sample(sample, args.fwithdrawn)
    logger.info('afr withdrawn {}'.format(len(sample)))

    if ethnic == 'eur':
        field_exclude = ['22019',  # chrom aneploidy
                         '22027'  # heterozygosity missing outlier
                         ]
        field_include = ['22006',  # white british
                         '22020',  # in pca = unrelated
                         ]
        field_exclude_val1 = []
    elif ethnic == 'noneur':
        # not white british
        field_exclude = ['22019',  # chrom aneploidy
                         '22027',  # heterozygosity missing outlier
                         ]
        field_exclude_val1 = [
            '22006',  # white british
        ]
        field_include = [
            '22020',  # in pca = unrelated
        ]
    else:
        raise RuntimeError

    fields_load = field_exclude + field_exclude_val1 + field_include
    # fields_load = field_exclude + field_include

    logger.debug('load ukb_whole')
    ukb = ukbfield.load_ukb(args.fukb_whole, fields_load)
    # ukb = ukb.rename(columns={'eid': 'id'})
    ukb = ukb.set_index('id')

    logger.debug('ukb\n{}'.format(ukb.head()))

    # make same order as qc_case.py => why??
    # 22006->22019->22027->22006
    # field = '22006'
    # logger.debug('include field {}'.format(field))
    # sample = include_or_exclude_sample_ukb_val1(sample, ukb, field)
    # logger.info('afr {}'.format(len(sample)))
    #
    # for field in field_exclude:
    #    logger.debug('exclude field {}'.format(field))
    #    sample = exclude_sample_ukb(sample, ukb, field)
    #    logger.info('afr {}'.format(len(sample)))
    #
    # field = '22020'
    # logger.debug('include field {}'.format(field))
    # sample = include_or_exclude_sample_ukb_val1(sample, ukb, field)
    # logger.info('afr {}'.format(len(sample)))

    for field in field_exclude:
        logger.debug('exclude field {}'.format(field))
        sample = exclude_sample_ukb(sample, ukb, field)
        logger.info('afr {}'.format(len(sample)))

    for field in field_exclude_val1:
        logger.debug('exclude_val1 field {}'.format(field))
        sample = include_or_exclude_sample_ukb_val1(sample, ukb, field, op='exclude')
        # sample = include_sample_ukb_val1(sample, ukb, field)
        # sample = include_or_exclude_sample_ukb_val1(sample, ukb, field)
        logger.info('afr {}'.format(len(sample)))

    for field in field_include:
        logger.debug('include field {}'.format(field))
        sample = include_or_exclude_sample_ukb_val1(sample, ukb, field)
        # sample = include_sample_ukb_val1(sample, ukb, field)
        logger.info('afr {}'.format(len(sample)))

    logger.debug('sample {}'.format(sample.head()))

    logger.info('write to {}'.format(args.fout))
    sample.index.to_series().to_csv(args.fout, header=False, index=False)


def extract_valid_id(sample):
    # sample = sample.loc[sample.index >= 0, :]
    sample = sample.loc[~sample.index.str.startswith('-'), :]
    return sample


def exclude_sample_ukb(sample, ukb, field):
    col_use = [col for col in ukb.columns if col.startswith(field)]
    if len(col_use) >= 2:
        raise RuntimeError('assume only one col')
    col_use = col_use[0]

    logger.debug('uniq {}'.format(ukb[col_use].unique()))

    # print('index', ukb[col_use].isnull())
    # exclude 1 = include not null
    ukb_include_index = ukb.loc[ukb[col_use].isnull()].index
    # print('index', ukb_include_index)

    assert sample.index.name == 'id'
    assert ukb.index.name == 'id'
    index_use = sample.index.join(ukb_include_index, how='inner')
    sample = sample.loc[index_use]
    return sample


def include_or_exclude_sample_ukb_val1(sample, ukb, field, op='include'):
    col_use = [col for col in ukb.columns if col.startswith(field)]
    if len(col_use) >= 2:
        raise RuntimeError('assume only one col')
    col_use = col_use[0]

    logger.debug('uniq {}'.format(ukb[col_use].unique()))

    if op == 'include':
        # include 1
        ukb_include_index = ukb.loc[ukb[col_use] == '1'].index
    elif op == 'exclude':
        # exclude 1
        ukb_include_index = ukb.loc[ukb[col_use] != '1'].index
    else:
        RuntimeError

    assert sample.index.name == 'id'
    assert ukb.index.name == 'id'
    index_use = sample.index.join(ukb_include_index, how='inner')
    sample = sample.loc[index_use]
    return sample


def exclude_sample(sample, fqc):
    cols = 'eid'.split(',')
    qc = pd.read_csv(fqc, sep=",",
                     names=cols,
                     header=None,
                     dtype=str)
    qc = qc.set_index('eid')

    index_intersection = sample.index.intersection(qc.index)

    sample = sample.drop(index=index_intersection)

    return sample


# def include_sample(sample, fqc):
#    cols = 'eid,qc'.split(',')
#    qc = pd.read_csv(fqc, sep=",",
#                     names=cols,
#                     header=0)
#    qc = qc.set_index('eid')
#    qc = qc.sort_index()
#
#    qc = qc[qc['qc'] == 1]
#
#    # here, index name of fam vanished
#    index_intersection = sample.index.intersection(qc.index)
#    sample = sample.loc[index_intersection]
#
#    return sample


def main():
    args = argument()

    if args.qc_sample:
        # logger.info('qc_sample')
        write_qc_sample(args)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
