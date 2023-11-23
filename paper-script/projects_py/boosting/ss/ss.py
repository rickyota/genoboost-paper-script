
import pandas as pd
from logging import getLogger

from ..io import filename as iof
from ...genetics.ss import ss_io as ssio

logger = getLogger(__name__)


def gmodels_ss():
    gmodels = ['add', 'dom', 'rec', 'hetonly']
    return gmodels


def load_ssd_raw(dresult, phe, cvi, gmodels=['add', 'dom', 'rec', 'hetonly'], mid_path=None):
    # gwas
    ssd = {}
    # gmodels = ['add', 'dom', 'rec', 'hetonly']
    for gmodel in gmodels:
        fss = iof.file_ss(dresult, phe, cvi, gmodel, mid_path)
        ss = ssio.load_ss(fss)
        # ss = extract_ss_paper(ss)
        ssd[gmodel] = ss
    return ssd


def merge_wgt_ssd(wgt_ss, ssd, how='inner', use_sid_tmp=False):
    # Does not assume ss-only variants exist
    # TMP
    ss_cols_use = ['id', 'sid', '-log10(P)', 'effect', 'p']
    # ss_cols_use = ['id', '-log10(P)', 'effect', 'p']
    for gmodel, ss in ssd.items():
        ss = ss.loc[:, ss_cols_use].copy()
        logger.debug('ss {}'.format(ss.head()))

        if gmodel == 'hetonly':
            for gmodel_merge in ['overdom', 'overrec']:
                wgt_ss = _merge_wgt_ss(wgt_ss, ss, gmodel_merge, how, use_sid_tmp)
                # wgt_ss = merge_wgt_ss(wgt_ss, ss, gmodel_merge, use_sid_tmp)

            # ss.columns = [x + '_overdom' for x in ss_cols_use]
            # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_overdom', how='inner')
            # ss.columns = [x + '_overrec' for x in ss_cols_use]
            # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_overrec', how='inner')
        else:
            wgt_ss = _merge_wgt_ss(wgt_ss, ss, gmodel, how, use_sid_tmp)

            # [ss.columns = [x + '_' + gmodel for x in ss_cols_use]
            # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how='inner')
        logger.debug('merge wgt_ss {}'.format(wgt_ss.head()))
    return wgt_ss


def _merge_wgt_ss(wgt_ss, ss, gmodel, how='inner', use_sid_tmp=False):

    if use_sid_tmp:
        ss = ss.drop(columns=['id'])
        ss = ss.rename(columns=lambda x: x if x == 'sid' else x + '_' + gmodel)
        wgt_ss = pd.merge(wgt_ss, ss, on='sid', how=how)
    else:
        # avoid duplicated 'sid'
        ss = ss.drop(columns=['sid'])
        ss = ss.rename(columns=lambda x: x if x == 'id' else x + '_' + gmodel)
        wgt_ss = pd.merge(wgt_ss, ss, on='id', how=how)

    # logger.debug('ss {}'.format(ss.head()))
    # logger.debug('ss {}'.format(ss.columns))

    # avoid contamination of ss.columns
    # ss = ss.copy()
    # ss.columns = [x + '_' + gmodel for x in ss.columns]
    # logger.debug('ss {}'.format(ss.head()))
    # if use_sid_tmp:
    #    wgt_ss = pd.merge(wgt_ss, ss, left_on='sid', right_on='sid_' + gmodel, how=how)
    #    # wgt_ss = pd.merge(wgt_ss, ss, left_on='sid', right_on='sid_' + gmodel, how='inner')
    # else:
    #    wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how=how)
    #    # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how='inner')

    # wgt_ss = wgt_ss.drop(columns=['id_' + gmodel, 'sid_' + gmodel])
    # wgt_ss=wgt_ss.drop(columns=['id_'+gmodel])

    return wgt_ss


def judge_genetic_model_gwas(ss):

    # ss=ss.copy()

    # col: '-log10(P)_add', dom, rec, (het)
    # to distinguish overdom and overrec, effect_overdom is necessary
    col_add = '-log10(P)_add'
    col_dom = '-log10(P)_dom'
    col_rec = '-log10(P)_rec'
    col_overdom = '-log10(P)_overdom'
    col_eff_overdom = 'effect_overdom'
    col_model = 'model_est_gwas'

    logger.debug('ss {}'.format(ss.head()))
    ss.loc[:, col_model] = 'NaN'
    if '-log10(P)_overdom' not in ss.columns:
        raise RuntimeError('use het')
        filt_add = (ss[col_add] > ss[col_dom]) & (ss[col_add] > ss[col_rec])
        ss.loc[filt_add, col_model] = 'add'
        filt_dom = (~filt_add) & (ss[col_dom] > ss[col_rec])
        ss.loc[filt_dom, col_model] = 'dom'
        filt_rec = (~filt_add) & (~filt_dom)
        ss.loc[filt_rec, col_model] = 'rec'
        assert (ss[col_model] != 'NaN').any()
    else:
        # do not forget to judge overdom or overrec
        # assume p of overdom = overrec
        # assume effect of overdom = overrec
        filt_add = (ss[col_add] > ss[col_dom]) & (ss[col_add] > ss[col_rec]) & (ss[col_add] > ss[col_overdom])
        assert (ss.loc[filt_add, col_model] == 'NaN').all()
        ss.loc[filt_add, col_model] = 'add'

        filt_dom = (~filt_add) & (ss[col_dom] > ss[col_rec]) & (ss[col_dom] > ss[col_overdom])
        assert (ss.loc[filt_dom, col_model] == 'NaN').all()
        ss.loc[filt_dom, col_model] = 'dom'

        filt_rec = (~filt_add) & (~filt_dom) & (ss[col_rec] > ss[col_overdom])
        assert (ss.loc[filt_rec, col_model] == 'NaN').all()
        ss.loc[filt_rec, col_model] = 'rec'

        # filt_overdom = (~filt_add) & (~filt_dom) & (~filt_rec) & (ss['effect_overdom'] > 0)
        filt_overdom = (~filt_add) & (~filt_dom) & (~filt_rec) & (ss[col_eff_overdom] >= 0)
        assert (ss.loc[filt_overdom, col_model] == 'NaN').all()
        ss.loc[filt_overdom, col_model] = 'overdom'

        filt_overrec = (~filt_add) & (~filt_dom) & (~filt_rec) & (ss[col_eff_overdom] < 0)
        # filt_overrec = (~filt_add) & (~filt_dom) & (~filt_rec) & (~filt_overdom)
        assert (ss.loc[filt_overrec, col_model] == 'NaN').all()
        ss.loc[filt_overrec, col_model] = 'overrec'

        # if one of -log10 is np.nan, this might happen
        assert (ss[col_model] != 'NaN').any()

    return ss
