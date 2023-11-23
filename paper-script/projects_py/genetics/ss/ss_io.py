
import os
from decimal import Decimal
import numpy as np
import pandas as pd
from ...genetics.snv import snv as snvop
from logging import getLogger

logger = getLogger(__name__)


def load_ss(fss, align_a1_alt=False, maf_thre=None):
    '''
    id: rsid
    sid: 'chrom':'pos'
    sida: 'chrom':'pos':'a1':'a2'
    '''
    # TMP; for condrec, condadd of cv1-
    if not os.path.exists(fss):
        logger.info("skip since fss does not exist: {}".format(fss))
        return None

    ss = pd.read_csv(fss, sep='\t', dtype={'P': str})
    ss = ss.rename(columns=str.lower)

    if align_a1_alt:
        ss = align_a1_to_alt_ss(ss)

    ss = snvop.add_sids(ss)

    ss.loc[:, 'effect'] = np.log(ss['or'])

    # ss = ss.rename(columns={'a1': 'A1', 'a2': 'A2'})
    # TMP
    # if 'dataset.12' in str(fss):
    # ss.loc[:, 'id'] = ss['sida']

    # def str_from_log10_decimal(x):
    #    # https://stackoverflow.com/questions/36092203/why-is-1e400-not-an-int
    #    dec = Decimal(10)**(Decimal(-1) * Decimal(x))
    #    if dec < 10**(-4):
    #        num_str = '{:.6e}'.format(dec)
    #    else:
    #        num_str = '{:.6f}'.format(dec)
    #    return num_str

    def p_decimal_from_str(x):
        # https://stackoverflow.com/questions/36092203/why-is-1e400-not-an-int
        return Decimal(x)

    def log10_p_from_str(x):
        return float(Decimal(x).log10())

    # error
    # def p_f128_from_str(x):
    #    return np.float128(Decimal(x))

    ss.loc[:, 'p_str'] = ss['p']

    # float
    # np.float128 prints out 0.0 but actually 1e-334. confirmed by np.log10
    ss.loc[:, 'p'] = ss['p_str'].astype(np.float128)
    ss.loc[:, '-log10(P)'] = -ss['p'].map(np.log10).astype(np.float64)

    if (ss['p'] == 0.0).any():
        logger.error('log10(P): {}'.format((ss['p'] == 0.0).any()))
        logger.error('log10(P): {}'.format(ss.loc[ss['p'] == 0.0, 'p']))
        raise RuntimeError

    # check if p and log10(p) work fine
    # if -log10(p) and decimal match, it's fine
    ss_check = ss[['p_str', 'p', '-log10(P)']].copy()
    ss_check.loc[:, 'p_f64'] = ss_check['p_str'].astype(np.float64)
    # decimal from str: actually you can run log from decimal 'p'
    ss_check.loc[:, '-log10(P)_decimal'] = -ss_check['p_str'].apply(lambda x: log10_p_from_str(x))
    # decimal type
    ss_check.loc[:, 'p_decimal'] = ss_check['p_str'].apply(lambda x: p_decimal_from_str(x))
    assert np.allclose(ss_check['-log10(P)'], ss_check['-log10(P)_decimal'])
    if (ss_check['p_f64'] == 0.0).any():
        logger.info('p_f64: {}'.format(ss_check.loc[ss_check['p_f64'] == 0.0, ['p', 'p_f64', 'p_decimal', '-log10(P)', '-log10(P)_decimal']]))

    if maf_thre is not None:
        raise NotImplementedError
        # raise RuntimeError('ny for dataset<=11; not necessary for dataset12')
        # TODO: use ukb_imp_com.frq for maf
        # print('limit maf', maf_thre)
        # fco_phe = fco.replace('%', phenotype)
        # gt = pd.read_hdf(fco_phe + '.obs.gt', key='gt')
        # ss = pd.merge(ss, gt[['vid', 'maf']], left_on='sida', right_on='vid', how='left')
        # ss = ss[ss['maf'] >= maf_thre].copy()

    logger.debug('ss len: {}'.format(len(ss)))

    return ss


def align_a1_to_alt_ss(ss):
    filt_unmatch = (ss['alt'] != ss['a1'])

    ss.loc[filt_unmatch, 'a1_tmp'] = ss.loc[filt_unmatch, 'a2']
    ss.loc[filt_unmatch, 'a2'] = ss.loc[filt_unmatch, 'a1']
    ss.loc[filt_unmatch, 'a1'] = ss.loc[filt_unmatch, 'a1_tmp']

    assert (ss['alt'] == ss['a1']).all()
    assert (ss['ref'] == ss['a2']).all()

    ss.loc[filt_unmatch, 'effect'] = -ss.loc[filt_unmatch, 'effect']
    if 'or' in ss.columns:
        ss.loc[filt_unmatch, 'or'] = 1.0 / ss.loc[filt_unmatch, 'or']
    ss.loc[filt_unmatch, 'l95'] = 1.0 / ss.loc[filt_unmatch, 'u95']

    return ss
