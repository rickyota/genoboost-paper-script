

import pandas as pd
import numpy as np
from logging import getLogger

from ..io import filename as iof
from ...genetics.snv import snv as snvop

logger = getLogger(__name__)


def load_wgt(fwgt):
    return pd.read_csv(fwgt, delim_whitespace=True, dtype=str)


# def load_nsnv_wgt(dout, method, phe, cvi, parad, kind, model):
def load_nsnv_wgt(dout, method, phe, cvi, parad, mid_path):
    fwgt = iof.file_wgt(dout, method, phe, cvi, parad, mid_path)
    # logger.debug('fwgt: {}'.format(fwgt))
    wgt = load_wgt(fwgt)

    if method == 'boosting':
        wgt = wgt[wgt['kind'] == 'SNV']
        wgt = wgt.drop_duplicates('var')
        nsnv = len(wgt)
    else:
        nsnv = len(wgt)
    return nsnv


def load_wgt_method(fwgt, method):
    if method == 'boosting':
        return load_wgt_boost(fwgt)
    else:
        return load_wgt_prev(fwgt, method)


def load_wgt_boost(fwgt):
    wgt = pd.read_csv(fwgt, delim_whitespace=True)
    # wgt = load_wgt(fwgt)
    wgt = wgt.fillna({'chrom': -1, 'pos': -1})
    wgt = wgt.astype({'chrom': np.int32, 'pos': np.int32})

    # wgt = wgt.rename(columns={'a1': 'A1', 'a2': 'A2'})

    wgt = snvop.add_sids(wgt)

    # rust
    wgt.loc[:, 'kind'] = wgt['kind'].str.lower()
    # wgt.loc[:, 'kind'] = wgt['kind'].replace({'cov': 'covariant'})
    wgt.loc[:, 'model'] = wgt['model'].str.lower()
    wgt.loc[:, 'model'] = wgt['model'].replace({'binary': 'cat2', 'free': 'cat3'})

    # logger.debug('wgt {}'.format(wgt.columns))

    # TODO: stay 'alpha', 'const' for cov, but unnecessary
    # if 'scorem' in wgt.columns:
    # allow logitnomissing
    if 'score2' in wgt.columns:
        # logitnomissing
        # remove alpha, const
        # filt = wgt['score0'].isna()
        # wgt.loc[filt, 'score0'] = wgt.loc[filt, 'alpha']
        wgt = wgt.drop(columns=['alpha', 'const'])
    else:
        # logitadd
        # cp alpha to effect for boosting_add
        wgt.loc[:, 'effect'] = wgt['alpha']

    if 'eps' not in wgt.columns:
        wgt.loc[:, 'eps'] = np.nan

    if 'eff_eps' not in wgt.columns:
        wgt.loc[:, 'eff_eps'] = np.nan

    # logger.debug('wgt {}'.format(wgt.columns))
    # logger.debug('wgt {}'.format(wgt.columns))

    # if 'alpha' in wgt.columns:
        # wgt.loc[:, 'effect'] = wgt['alpha']
    # logger.debug('wgt {}'.format(wgt.columns))

    return wgt


# TOFIX: should be done in format_wgt
def load_wgt_prev(fwgt, method):
    if method == 'ldpred':
        # cols = 'chrom,pos,sid,A1,A2,N1,effect'.split(',')
        # cols = 'chrom,pos,sid,A1,A2,N1,effect,N2'.split(',')
        # wgtldpred = pd.read_csv(fwgt_ldpred, sep='\\s+', header=0,
        # wgtldpred = pd.read_csv(fwgt_ldpred, sep='\t', header=0,
        #                        names=cols)
        wgt = pd.read_csv(fwgt, sep='\t', header=0)
        if len(wgt.columns) == 7:
            # wgt.columns = 'chrom,pos,sid,a1,a2,N1,effect'.split(',')
            wgt.columns = 'chrom,pos,id,a1,a2,N1,effect'.split(',')
        elif len(wgt.columns) == 8:
            # wgt.columns = 'chrom,pos,sid,a1,a2,N1,effect,N2'.split(',')
            wgt.columns = 'chrom,pos,id,a1,a2,N1,effect,N2'.split(',')
        else:
            raise RuntimeError("ldpred cols length was ", len(wgt.columns))
        # print("wgtldrped", wgtldpred)
        # print("wgtldrped", wgtldpred.info())
        # print("wgtldrped[A2]", wgtldpred['A2'])
    elif method == 'clump':
        # cols = 'chrom,pos,sid,A1,A2,effect'.split(',')
        # cols = 'chrom,pos,sid,a1,a2,effect,p'.split(',')
        cols = 'chrom,pos,id,a1,a2,effect,p'.split(',')
        wgt = pd.read_csv(fwgt, delim_whitespace=True, header=0, names=cols)
    elif method == 'snpnet':
        # cols = 'sida,chrom,pos,a1,a2,effect'.split(',')
        cols = 'id,chrom,pos,a1,a2,effect'.split(',')
        wgt = pd.read_csv(fwgt, delim_whitespace=True, header=0, names=cols)
    elif method == 'lassosum':
        # cols = 'chrom,sid,pos,a1,a2,effect'.split(',')
        cols = 'chrom,id,pos,a1,a2,effect'.split(',')
        wgt = pd.read_csv(fwgt, delim_whitespace=True, header=0, names=cols)
    elif method == 'sbayesr':
        # cols = 'rs,chrom,pos,a1,a2,effect'.split(',')
        cols = 'id,chrom,pos,a1,a2,effect'.split(',')
        wgt = pd.read_csv(fwgt, delim_whitespace=True, header=0, names=cols)
    elif method == 'prscs':
        # cols = 'chrom,rs,pos,a1,a2,effect'.split(',')
        cols = 'chrom,id,pos,a1,a2,effect'.split(',')
        wgt = pd.read_csv(fwgt, delim_whitespace=True, header=0, names=cols)
    elif method == 'snpboost':
        cols = 'id,chrom,pos,a1,a2,effect'.split(',')
        wgt = pd.read_csv(fwgt, delim_whitespace=True, header=0, names=cols)
    else:
        raise NotImplementedError('Unknown method', method)

    wgt = snvop.add_sids(wgt)

    return wgt


def estimate_model_score_boost(wgt_snv_find, score1='score1_rel', score2='score2_rel', genetic_threshold_kind='1'):
    # should include 'score1' and 'score2'
    # this is relative score to 'score0'

    wgt_snv_find.loc[:, 'ratio'] = wgt_snv_find[score1] / wgt_snv_find[score2]

    # print('ratio nan',wgt_snv_find.loc[wgt_snv_find['ratio'].isna(),'ratio'].to_string())

    # for 'ratio'==np.nan
    # -> cov
    if wgt_snv_find['ratio'].isnull().any():
        print(wgt_snv_find.loc[wgt_snv_find['ratio'].isnull(), :])
        print(wgt_snv_find.loc[wgt_snv_find['ratio'].isnull(), ['ratio', score1, score2]])
        raise RuntimeError('ratio=nan is in wgt')

    # bins = [-np.inf, -0.249, 0.251, 0.75, 1.249, np.inf]
    # bins = [-np.inf, -0.249, 0.249, 0.75, 1.249, np.inf]

    if genetic_threshold_kind == '1':
        bins = [-np.inf, -0.25, 0.25, 0.75, 1.25, np.inf]
    elif genetic_threshold_kind == '2':
        # 1/2.828, 1/1.414
        bins = [-np.inf, -0.354, 0.353, 0.708, 1.414, np.inf]
    else:
        raise NotImplementedError('Unknown genetic_threshold_kind: ', genetic_threshold_kind)

    labels = ['hetonly', 'rec', 'add', 'dom', 'hetonly']
    # ordered=False for duplicated labels
    wgt_snv_find.loc[:, 'model_est'] = pd.cut(wgt_snv_find['ratio'], bins=bins, labels=labels,
                                              ordered=False, right=True, include_lowest=True).astype('string')
    # wgt_snv_find.loc[:, 'model_est'] = pd.cut(wgt_snv_find['ratio'], bins=bins, labels=labels, ordered=False).astype('string')
    # logger.debug('model_est {}'.format(wgt_snv_find.to_string()))

    # if hetonly and rs1 > 0; overdom
    def fn_hetonly(x):
        # logger.debug('x {}'.format(x))
        if x['model_est'] != 'hetonly':
            return x['model_est']

        if x[score1] >= 0:
            return 'overdom'
        else:
            return 'overrec'

    wgt_snv_find.loc[:, 'model_est'] = wgt_snv_find[[score1, 'model_est']].apply(fn_hetonly, axis=1)

    # bins = [-np.inf, -0.25, 0.25, 0.75, 1.25, np.inf]
    # labels = ['overrec', 'rec', 'add', 'dom', 'overdom']
    # wgt_snv_find.loc[:, 'model_est'] = pd.cut(wgt_snv_find['ratio'], bins=bins, labels=labels).astype('string')

    return wgt_snv_find


def format_wgt_boost(wgt, nsnv=None, genetic_threshold_kind='1'):
    # TODO: rename fn
    '''
    Dedup snv.
    Extract first nsnv.
    '''

    wgt = wgt[wgt['kind'] == 'snv'].copy()

    # wgt.loc[:, 'sid'] = wgt['boost_sid']
    # wgt.loc[:, 'vid'] = wgt['var']
    wgt = wgt.rename(columns={'var': 'id'})

    # calculate relative score

    def add_relative_score(wgt):
        wgt.loc[:, 'score1_rel'] = wgt['score1'] - wgt['score0']
        wgt.loc[:, 'score2_rel'] = wgt['score2'] - wgt['score0']
        return wgt

    if 'score0' in wgt:
        wgt = add_relative_score(wgt)
        # add rec or dom
        wgt = estimate_model_score_boost(wgt, score1='score1_rel', score2='score2_rel', genetic_threshold_kind=genetic_threshold_kind)

        # in paper.py
        # wgt.loc[:, 'alpha'] = wgt['score1_rel']
        # wgt.loc[:, 'alpha'] = wgt['score2_rel']
        # TMP: make score2=alpha
        # wgt.loc[:, 'alpha'] = wgt['score2']

        wgt.loc[:, 'dom'] = wgt['score1'] - wgt['score0']
        wgt.loc[:, 'rec'] = wgt['score2'] - wgt['score1']

        # add assuming 2 * score1_rel = score2_rel
        # score1 - mean(score0, score2)
        wgt.loc[:, 'add'] = (wgt['dom'] + wgt['rec']) / 2
        # overdom assuming score0~score2
        # score1 - mean(score0, score2)
        wgt.loc[:, 'overdom'] = wgt['score1_rel'] - (0.0 + wgt['score2_rel']) / 2
        # overrec is the same
        wgt.loc[:, 'overrec'] = wgt['overdom']

    cols_unchange = 'id,sid,sida,chrom,pos,a1,a2,effect'.split(',')
    wgt = wgt.rename(columns=lambda x: 'boost_' + x if x not in cols_unchange else x)
    # wgt.columns = ['boost_' + col for col in wgt.columns]
    # wgt = wgt.rename(columns={'boost_id': 'id',
    #                          'boost_sid': 'sid',
    #                          'boost_sida': 'sida',
    #                          'boost_chrom': 'chrom',
    #                          'boost_pos': 'pos',
    #                          'boost_a1': 'a1',
    #                          'boost_a2': 'a2',
    #                          'boost_effect': 'effect',  # for boosting_add
    #                          })

    # since niter_boost is for unduplicated snvs
    wgt = wgt.drop_duplicates(subset=['id'], keep='first')

    # wgt = wgt.drop_duplicates(subset=['boost_id'], keep='first')

    if nsnv is not None:
        wgt = wgt.iloc[:nsnv, :]

        # print("wgt after niter_pred")
        # print("len(wgt)", len(wgt))
        # print("wgt", wgt)

    return wgt
