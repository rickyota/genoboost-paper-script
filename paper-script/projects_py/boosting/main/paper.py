# model = key of method_model
# ex. 'boosting_nonadd': method='boosting', model='logitnomissing_...'
# ex. 'snpnet_nonadd': method='snpnet', model='nonadd'

import argparse
import os
import pathlib
import yaml
import functools
# from decimal import Decimal
from collections import defaultdict
import numpy as np
import pandas as pd
from sklearn import metrics
from scipy import stats as scistats
from logging import getLogger

from ...system.logger import logger_setting
from ...genetics.sample import sampleio as sampleio
from ...genetics.sample import sample as sampleop
from ...genetics.score import scoreio as scoreio
from ...genetics.snv import snv as snvop
# from ...genetics.ss import ss_io as ssio
from ..io import filename as iof
from ..io import arg as ioarg
# from ..score import score as scoreop
from ..wgt import wgt as wgtop
from ..ss import ss as ssop
from ..paper import validate as paper_validate
# for paper public
#from ..paper import config as paper_config
from ..paper import public_stat as paper_public_stat

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    # FIXME: add --boost-no-monitor
    # since with-monitor is default

    parser.add_argument('--runs',
                        nargs="*",
                        required=True,
                        help='')

    parser.add_argument('--format',
                        nargs="*",
                        help='')

    parser.add_argument('--dout',
                        help='output dir')

    # ./result/nonaddpgs/dataset.12/
    parser.add_argument('--dresult',
                        help='result dir')

    parser.add_argument('--ddata',
                        help='data dir')

    # assume  --method-mid-path $(method) $(model)
    parser.add_argument('--method-mid-path',
                        nargs="+",
                        action='append',
                        help='')

    # assume  --method-mode $(method) $(model)
    parser.add_argument('--method-mode',
                        nargs="+",
                        action='append',
                        help='')

    # assume  --model-integrate boosting_integrate "boosting_add boosting_nonadd"
    parser.add_argument('--model-integrate',
                        nargs="+",
                        action='append',
                        help='')

    parser.add_argument('--methods',
                        nargs="*",
                        # help='Do not use "_" for method',
                        )

    parser.add_argument('--methods-type',
                        nargs="*",
                        help='')

    # TODO: --methods-type "noneur" "boosting_integrate_noneur snpnet_noneur"

    # use --para-regon
    # parser.add_argument('--regon-ss',
    #                    help='')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--fcv',
                        help='')

    parser.add_argument('--cvn',
                        help='')

    parser.add_argument('--phes',
                        nargs="*",
                        help='')

    parser.add_argument('--phes-female',
                        nargs="*",
                        help='')

    parser.add_argument('--para-cand',
                        help='')

    parser.add_argument('--para-com',
                        help='')

    parser.add_argument('--para-regon',
                        help='')

    parser.add_argument('--fgene-exon',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def get_fymld(ddata, para_cand, para_regon, para_com):
    dpara = ddata / 'para_paper'
    fregon = dpara / ('para.' + para_regon + '.yml')
    # fregon_nocov = './para.' + para_regon + '.nocov.yml'
    fcand = dpara / ('para.' + para_cand + '.yml')
    fcom = dpara / ('para.' + para_com + '.yml')
    # fcond = './para.' + para_cond + '.yml'

    fymld = {
        'regon': fregon,
        # 'regon_nocov': fregon_nocov,
        'cand': fcand,
        'com': fcom,
        # 'condition': fcond
    }
    return fymld


# when withcov=False, regon is not necessary
# def load_allstat(dresult, method, phe, withcov, regon=None, kind=None, model=None):
@functools.cache
def load_allstat(dresult, method, phe, withcov, regon=None, mid_path=None, mode=None,):

    # (kind, model) = mid_path
    fallstat = iof.file_allstat(dresult, method, phe, withcov, regon, mid_path, mode)
    logger.debug('fallstat: {}'.format(fallstat))
    # print('fallstat',fallstat)
    allstat = pd.read_csv(fallstat, sep='\t', index_col=0)
    # FIXME: para should be string, otherwise, cannot access to filename
    # -> but troublesome for cvi, nsnv etc.
    # allstat = pd.read_csv(fallstat, sep='\t', index_col=0, dtype=str)

    # print('allstat',allstat)
    return allstat


def load_allstatd_integrate(dresult, method, phe, withcov, method_modeld, methodd_mode, models_integrate, regon=None):
    # methodm='boosting_integrate'
    assert method == 'boosting'
    # assert methodm.startswith('boosting_integrate')

    allstatd = {}
    # for method_genet in model_integrate[methodm]:
    # for method_genet in ['boosting_add', 'boosting_nonadd']:
    for method_genet in models_integrate:
        if methodd_mode is None:
            mode = None
        else:
            mode = methodd_mode[method_genet]
        allstat = load_allstat(dresult, method, phe, withcov, regon, method_modeld[method_genet], mode)
        if len(allstat) == 0:
            logger.warning('alsltat len=0: {}'.format(method))
        # print('method_genot',method_genet)
        # print('allstat,',allstat)
        allstatd[method_genet] = allstat

    return allstatd


def convert_string_none(d):
    # convert 'None' to None
    # max depth=3 now

    for k, v in d.items():
        if v == 'None':
            d[k] = None
        if isinstance(v, dict):
            for k2, v2 in v.items():
                if v2 == 'None':
                    d[k][k2] = None
                if isinstance(v2, dict):
                    for k3, v3 in v2.items():
                        if v3 == 'None':
                            d[k][k2][k3] = None
    return d


def load_para_yml_classic(f_para_yml):
    """
    load yml config
    """

    with open(f_para_yml) as file:
        yml = yaml.safe_load(file)
        # load as string
        # yml = yaml.safe_load(file, Loader=yaml.BaseLoader)
    return yml


@functools.cache
def load_regon_yml(fyml):
    yml = load_para_yml_classic(fyml)
    regon = convert_string_none(yml)
    # print('load_regon_yml', regon)
    return regon


def load_regon(fymld):
    fregon = fymld['regon']
    regond = load_regon_yml(fregon)
    return regond


@functools.cache
def load_cand(fyml):
    yml = load_para_yml_classic(fyml)
    yml = convert_string_none(yml)
    # print('load_cand_yml', yml)
    nsnvs = yml['nsnvs']
    # if 'lr' in yml:
    lrs = yml['lr']
    # else:
    #    # do not use lr
    #    # TODO: want to change
    #    lrs = ()
    return nsnvs, lrs


def extract_cand_allstat(allstat, cvi, method, va_on, fcand):

    nsnvs, lrs = load_cand(fcand)

    # logger.debug('alltat in bestpara: {}'.format(allstat))
    allstat_debug = allstat.copy()

    allstat = allstat.loc[allstat['dataset'] == va_on, :]
    allstat = allstat.loc[allstat['cvi'] == cvi, :]

    # logger.debug('alltat in bestpara extract: {}'.format(allstat))

    if len(allstat) == 0:
        logger.debug('alltat in bestpara: {}'.format(allstat_debug))
        logger.warning('WARNING: dataset is not in allstat: {}, {}, {}'.format(method, cvi, va_on))
        return None

    # limit nsnvs
    if method == 'boosting':
        if nsnvs is not None:
            # -1 for monitor
            nsnvs = nsnvs + [-1]
            allstat = allstat.loc[allstat['n'].isin(nsnvs), :]
            # allstat = allstat.loc[allstat['n'].isin(nsnvs), :]
    elif method in iof.methods_classic_use_nsnvs():
        if nsnvs is not None:
            # allow all of more than largest
            max_nsnv = max(nsnvs)
            # FIXME: if allstat does not have cols defined in iof, add null somewhere
            # TMP if
            # if 'n' in allstat:
            allstat = allstat.loc[allstat['n'].isin(nsnvs) | (allstat['n'] >= max_nsnv) | (allstat['n'] == -1), :]
            # else:
            #    # TMP
            #    logger.warning('n is not in allstat')
            #    allstat.loc[:, 'n'] = -1

    if method == 'boosting':
        assert 'lr' in allstat.columns
        # None : use the same as in config
        if (lrs != ()) and (lrs is not None):
            lrs = lrs + [-1]
            # ng: this overwrite lrs
            # lrs += [-1]
            allstat = allstat.loc[allstat['lr'].isin(lrs), :]
            # allstat = allstat.loc[allstat['lr'].isin(lrs), :]
        # raise NotImplementedError
        # if 'boosting_learning_rate' in allstat.columns:
        #    # None : use the same as in config
        #    if (lrs != ()) and (lrs is not None):
        #        allstat = allstat.loc[allstat['boosting_learning_rate'].isin(lrs), :]

    # ok
    # print('lrs in allstat',allstat['lr'].unique())
    # print('allstat',allstat)

    # limit allstat for args
    # if n is not None:
    #    if smaller:
    #        allstat = allstat.loc[allstat['nsnvs'] <= n, :]
    #    else:
    #        allstat = allstat.loc[allstat['nsnvs'] == n, :]

    allstat = allstat.copy()

    return allstat


def load_best_para(allstat, cvi, method, max_on, va_on, fcand,
                   return_stat=False):

    allstat = extract_cand_allstat(allstat, cvi, method, va_on, fcand)

    if allstat is None:
        return None

    # nsnvs, lrs = load_cand(fcand)

    # print('alltat in bestpara',allstat)

    # allstat = allstat.loc[allstat['dataset'] == va_on, :]
    # allstat = allstat.loc[allstat['cvi'] == cvi, :]

    # if len(allstat) == 0:
    #    logger.warning('WARNING: dataset is not in allstat: {}, {}, {}'.format(method, cvi, va_on))
    #    return None

    # limit nsnvs
    # if method == 'boosting':
    #    if nsnvs is not None:
    #        allstat = allstat.loc[allstat['n'].isin(nsnvs), :]
    # elif method in iof.methods_classic_use_nsnvs():
    #    if nsnvs is not None:
    #        # allow all of more than largest
    #        max_nsnv = max(nsnvs)
    #        allstat = allstat.loc[allstat['n'].isin(nsnvs) | (allstat['n'] >= max_nsnv) | (allstat['n'] == -1), :]

    # if method == 'boosting':
    #    assert 'lr' in allstat.columns
    #    # None : use the same as in config
    #    if (lrs != ()) and (lrs is not None):
    #        allstat = allstat.loc[allstat['lr'].isin(lrs), :]
    #    # raise NotImplementedError
    #    # if 'boosting_learning_rate' in allstat.columns:
    #    #    # None : use the same as in config
    #    #    if (lrs != ()) and (lrs is not None):
    #    #        allstat = allstat.loc[allstat['boosting_learning_rate'].isin(lrs), :]

    # ok
    # print('lrs in allstat',allstat['lr'].unique())
    # print('allstat',allstat)

    # limit allstat for args
    # if n is not None:
    #    if smaller:
    #        allstat = allstat.loc[allstat['nsnvs'] <= n, :]
    #    else:
    #        allstat = allstat.loc[allstat['nsnvs'] == n, :]

    paras_method = iof.get_paras_method(method)
    cols_para = ['nsnv_use'] + paras_method
    # cols_para = ['n', 'nsnv_use'] + paras_method
    allstat = allstat.set_index(cols_para)
    stat_nsnv_idx_ = allstat[max_on].idxmax(skipna=True)
    allstat_max = allstat.loc[[stat_nsnv_idx_], :].reset_index(drop=False).iloc[0]
    stat_best = allstat_max.loc[max_on]

    # TODO: using str is the best, but troublesome in ex. allstst['n']>=max_nsnv
    paras_best_method = {}
    for para in paras_method:
        paras_best_method[para] = allstat_max.loc[para]
        # n=-1 should be fine? since later could be in new df
        # if n == -1:
        #    n=None
        # else:
        #    n=n
    # 'n' is in para for all allstat
    # paras_best_method['n'] = allstat_max.loc['n']
    paras_best_method['cvi'] = allstat_max.loc['cvi']
    # nsnvs_fname = allstat_max.loc['nsnvs_fname']
    # if nsnvs_fname==-1:
    #    nsnvs_fname = None
    # else:
    #    nsnvs_fname = nsnvs_fname
    paras_best_method['nsnv_use'] = allstat_max.loc['nsnv_use']

    if return_stat:
        return paras_best_method, stat_best
    else:
        return paras_best_method


def load_best_para_integrate(allstatd, cvi, method, models_integrate, max_on, va_on, fcand, return_stat=False):

    paras_best = None
    stat_best = None
    model_best = None

    for model in models_integrate:
        paras_best_in_model, stat_best_in_model = load_best_para(allstatd[model], cvi, method, max_on, va_on, fcand, return_stat=True)
        if stat_best is None or stat_best_in_model > stat_best:
            paras_best = paras_best_in_model
            stat_best = stat_best_in_model
            model_best = model

    assert stat_best is not None

    paras_best |= {'method_integrate': model_best}
    if return_stat:
        return paras_best, stat_best
    else:
        return paras_best

    # paras_best_add, stat_best_add = load_best_para(allstatd['boosting_add'], cvi, method, max_on, va_on, fcand, return_stat=True)
    # print('paras_best_add',paras_best_add)
    # paras_best_nonadd, stat_best_nonadd = load_best_para(allstatd['boosting_nonadd'], cvi, method, max_on, va_on, fcand, return_stat=True)

    # if stat_best_add > stat_best_nonadd:
    #    paras_best_add |= {'method_integrate': 'boosting_add'}
    #    if return_stat:
    #        return paras_best_add, stat_best_add
    #    else:
    #        return paras_best_add
    # else:
    #    paras_best_nonadd |= {'method_integrate': 'boosting_nonadd'}
    #    if return_stat:
    #        return paras_best_nonadd, stat_best_nonadd
    #    else:
    #        return paras_best_nonadd


def load_stat_from_para_best(allstat,
                             paras_best,
                             dataset,
                             stat_on,
                             method):

    allstat = allstat.copy()
    allstat = allstat.loc[allstat['dataset'] == dataset, :]

    paras_method = iof.get_paras_method(method)

    # print('paras_best', paras_best)

    # col_exclude = ['cvi', 'n']
    # col_exclude = ['cvi', 'n', 'nsnvs_fname']
    # col_index = ['cvi', 'nsnvs', 'nsnvs_fname']
    # col_index = ['cvi', 'nsnvs', 'nsnvs_fname']

    # paras_col = [col for col in paras_best if col not in col_exclude]
    # cols_para = col_index + [method + '_' + col for col in paras_col]

    cols_para = ['cvi'] + paras_method
    # cols_para = ['cvi', 'n', paras_method]

    # nsnvs_fname is not para here
    idx = tuple(paras_best[col] for col in cols_para)
    # idx = (paras_best['cvi'], paras_best['n']) + tuple(paras_best[col] for col in paras_method)
    # idx = (paras_best['cvi'], paras_best['n'], paras_best['nsnvs_fname']) + tuple(paras_best[col] for col in paras_col)

    # print('idx, cols_para', idx, cols_para)

    allstat = allstat.set_index(cols_para)

    if stat_on not in allstat:
        raise RuntimeError('stat not in allstat: {} in {}'.format(stat_on, method))

    # print('allstat', allstat)
    stat = allstat.loc[idx, stat_on]

    return stat


def load_best_stat_from_allstat(allstat, paras_best, dataset, stat_on, method, methodm, phe, max_on):
    if paras_best is None:
        return None

    statd = {}
    stat = load_stat_from_para_best(
        allstat,
        paras_best,
        dataset,
        stat_on,
        method)
    statd |= {
        # store 'boosting_nonadd' not 'boosting'
        'method': methodm,
        # 'method': method,
        'phe': phe,
        'dataset': dataset,
        'stat': stat,
        'stat_on': stat_on,
        'max_on': max_on,
        'cvi': paras_best['cvi'],
        # some method does not have 'n' in paras_best
        # 'n': paras_best['n'],
        'nsnv_use': paras_best['nsnv_use'],
    }

    return statd


def load_best_stat_from_allstat_integrate(allstatd, paras_best, dataset, stat_on, method, methodm, phe, max_on):
    if paras_best is None:
        return None

    method_integrate = paras_best['method_integrate']
    allstat = allstatd[method_integrate]

    statd = load_best_stat_from_allstat(allstat, paras_best, dataset, stat_on, method, methodm, phe, max_on)

    # TODO:
    # statd|={'method_integrate':method_integrate}

    return statd


def load_best_stat_from_allstat_cvs_datasets(
        dresult,
        methodm,
        phe,
        withcov,
        reg_on,
        stat_on,
        max_on,
        cvn,
        method_modeld,
        methodd_mode_ts,
        model_integrate,
        paras_best_cvs,
        datasets):
    # 'boosting_nonadd' -> 'boosting'
    method = methodm_to_method(methodm)

    stats = []
    # if methodm == 'boosting_integrate':
    if methodm.startswith('boosting_integrate'):
        # allstat for stat
        allstatd = load_allstatd_integrate(dresult, method, phe, withcov, method_modeld, methodd_mode_ts, model_integrate[methodm], reg_on)
        for cvi in range(cvn):
            paras_best = paras_best_cvs[cvi]
            for dataset in datasets:
                statd = load_best_stat_from_allstat_integrate(allstatd, paras_best, dataset, stat_on, method, methodm, phe, max_on)
                stats.append(statd)
    else:
        # None if not boosting
        # model = method_modeld[methodm]

        if methodd_mode_ts is None:
            mode_ts = None
        else:
            mode_ts = methodd_mode_ts[methodm]

        # allstat for stat
        # try:
        allstat = load_allstat(dresult, method, phe, withcov, reg_on, method_modeld[methodm], mode_ts)
        # except:
        # could not read file

        for cvi in range(cvn):
            paras_best = paras_best_cvs[cvi]
            for dataset in datasets:
                statd = load_best_stat_from_allstat(allstat, paras_best, dataset, stat_on, method, methodm, phe, max_on)
                stats.append(statd)
    return stats


def methodm_to_method(methodm):
    if '_' in methodm:
        splits = methodm.split('_')
        if len(splits) != 2:
            raise RuntimeError('Wrong methodm. "_" should be once in methodm: ', methodm)
        method = splits[0]
    else:
        method = methodm
    # if methodm.startswith('boosting_'):
    #    method = 'boosting'
    # else:
    #    method = methodm
    return method


def load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on, va_on, fcand,
                       method_modeld, methodd_mode_va, model_integrate, cvn=None, cvis=None):
    # 'boosting_nonadd' -> 'boosting'
    method = methodm_to_method(methodm)

    if cvn is not None:
        cvis = list(range(cvn))

    paras_best_cvs = {}
    # if methodm == 'boosting_integrate':
    if methodm.startswith('boosting_integrate'):
        # allstat for maxon
        allstatd_maxon = load_allstatd_integrate(dresult, method, phe, withcov_at_maxon, method_modeld,
                                                 methodd_mode_va, model_integrate[methodm], reg_on)
        # for cvi in range(cvn):
        for cvi in cvis:
            paras_best = load_best_para_integrate(allstatd_maxon, cvi, method, model_integrate[methodm], max_on, va_on, fcand)
            # print('paras_best',paras_best)
            # paras_best = load_best_para_integrate(allstatd_maxon, cvi, method, max_on, va_on, fymld['cand'])
            paras_best_cvs[cvi] = paras_best
    else:

        if methodd_mode_va is None:
            mode = None
        else:
            mode = methodd_mode_va[method]

        # allstat for maxon
        allstat_maxon = load_allstat(dresult, method, phe, withcov_at_maxon, reg_on, method_modeld[methodm], mode)
        # allstat_maxon = load_allstat(dresult, method, phe, withcov_at_maxon, reg_on, method_modeld[methodm], methodd_mode_va[method])
        # for cvi in range(cvn):
        for cvi in cvis:
            paras_best = load_best_para(allstat_maxon, cvi, method, max_on, va_on, fcand)
            # paras_best = load_best_para(allstat_maxon, cvi, method, max_on, va_on, fymld['cand'])
            paras_best_cvs[cvi] = paras_best
    return paras_best_cvs


def regon_vaon(methodm, regond, withcov_at_maxon):
    # 'boosting_nonadd' -> 'boosting'
    method = methodm_to_method(methodm)

    if withcov_at_maxon:
        reg_on = regond[method]['regon']
    else:
        reg_on = None
    _, va_on = iof.get_datasets_vaon(reg_on)
    return reg_on, va_on


# phes_female not necessary here
# or split into load_para / load_best
def load_best_stats_acc(
    dresult,
    methodms,
    method_modeld,
    methodd_mode_va, methodd_mode_ts,
    model_integrate,
    phes,
    cvn,
    fymld,
    max_on,
    stat_on,
    withcov_at_maxon,
    withcov,
    datasets,
        method_modeld_va=None):

    if method_modeld_va is None:
        method_modeld_va = method_modeld.copy()

    regond = load_regon(fymld)
    stats = []
    for methodm in methodms:
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)

        # TMP!!
        if methodm.endswith('noneur'):
            reg_on, va_on = 'va', 'va'

        for phe in phes:
            # for allstat warning
            # logger.debug('phe: {}'.format(phe))
            paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                                va_on, fymld['cand'], method_modeld_va, methodd_mode_va, model_integrate, cvn)
            stats_cvs_datasets = load_best_stat_from_allstat_cvs_datasets(
                dresult, methodm, phe, withcov, reg_on, stat_on, max_on, cvn,
                method_modeld, methodd_mode_ts,
                model_integrate, paras_best_cvs, datasets)

            stats += stats_cvs_datasets

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)
    return stats


def load_phes_simu(ddata, simu_type):
    # FIXME:
    fphe = ddata / "simu" / simu_type / "phe" / "phe.phe"
    df = pd.read_csv(fphe, sep='\t', header=0)
    phes = df.columns.tolist()
    phes.remove('IID')
    return phes


def add_stat_simu(stats, simu_type):
    # TODO: fn
    simus = simu_type.split('_')
    for stat in simus:
        if '-' not in stat:
            continue
        x = stat.split('-')
        stats[x[0]] = x[1]

    stats['h2'] = str(float(stats['h2add']) + float(stats['h2dom']))
    return stats


def load_best_stats_acc_simu(
    ddata,
        dresult,
        methodms,
        method_modeld,
        model_integrate,
        simu_types,
        # phes,
        cvn,
        fymld,
        max_on,
        stat_on,
        withcov_at_maxon,
        withcov,
        datasets):
    regond = load_regon(fymld)
    stats = []
    for methodm in methodms:
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
        for simu_type in simu_types:

            phes = load_phes_simu(ddata, simu_type)
            # add simu_type here
            method_modeld_simu = method_modeld.copy()
            # logger.debug('method_modeld_simu: {}'.format(method_modeld_simu))
            # methods_keys= list(method_modeld_simu.keys())
            # if methodm in methods_keys:
            for k in method_modeld_simu.keys():
                v = method_modeld_simu[k]
                v = list(v)
                # boosting: 'simu' / simu_type / model
                # ldrped: 'simu' / simu_type
                v.insert(1, simu_type)
                v = tuple(v)
                method_modeld_simu[k] = v
            # TODO: not to use defaultkey?
            # else:
            #    method_modeld_simu[methodm] =('default',simu_type)

            for phe in phes:
                # for allstat warning
                # logger.debug('phe: {}'.format(phe))
                paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                                    va_on, fymld['cand'], method_modeld_simu, None, model_integrate, cvn)
                stats_cvs_datasets = load_best_stat_from_allstat_cvs_datasets(
                    dresult, methodm, phe, withcov, reg_on, stat_on, max_on, cvn,
                    method_modeld_simu, None,
                    model_integrate, paras_best_cvs, datasets)

                stats_cvs_datasets = [add_stat_simu(x, simu_type) for x in stats_cvs_datasets]
                # stats_cvs_datasets = add_stat_simu(stats_cvs_datasets, simu_type)

                stats += stats_cvs_datasets

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)
    return stats


def load_score(fscore, paras_best, cv, pheno, phe, dataset, sex):
    score_whole = scoreio.load_score_pheno(fscore, pheno, phe)

    if 'n' not in paras_best or paras_best['n'] == -1:
        col = 'score'
    else:
        col = 'score_n-' + str(paras_best['n'])

    score_col = score_whole.loc[:, ['id', 'status', col]].rename(columns={col: 'score'}).copy()
    score = sampleop.extract_dataset(score_col, dataset, paras_best['cvi'], cv, sex)
    if len(score) == 0:
        raise RuntimeError('#samples=0.')
    return score


def calc_or(score, or_props):
    n = len(score)
    case = (score['status'] == 1).sum()
    cont = (score['status'] == 0).sum()
    sc_sorted = score.sort_values(by='score', ascending=False, ignore_index=True)

    statds = []
    for prop in or_props:
        high = int(n * prop * 0.01)
        thre = sc_sorted['score'].iloc[high]
        # cannot use score.iloc[:high] since many samples could have the same scores
        # should be '=' since > might result in casehigh=0
        highs = sc_sorted.loc[(sc_sorted['score'] >= thre), :]
        # print("prop vs actual prop", high, len(highs))
        if np.abs(high - len(highs)) > high * 0.2:
            logger.warning("warning: prop vs actual prop: {} {} in {}%".format(high, len(highs), prop))
            # logger.warning("warning: prop vs actual prop: {} {} in {} {} {}%".format(high, len(highs), method, phe, prop))
        casehigh = (highs['status'] == 1).sum()
        conthigh = (highs['status'] == 0).sum()
        caselow = case - casehigh
        contlow = cont - conthigh
        or_stat = (casehigh * contlow) / (conthigh * caselow)
        # print("or_",or_)

        statd = {
            'stat_on': 'or_prop-' + str(prop),
            'stat': or_stat
        }
        statds.append(statd)

    return statds


def calc_best_stat_or(score, methodm, paras_best, dataset, phe, max_on, or_props):
    if paras_best is None:
        return []

    # for warning
    logger.debug('method, phe: {}, {}'.format(methodm, phe))
    # get statds: or_prop3 etc.
    statds = calc_or(score, or_props)

    statd_para = {
        # TODO? 'methodm': method
        # but troublesome in plot
        'method': methodm,
        # 'method': method,
        'phe': phe,
        'dataset': dataset,
        # in statds
        # 'stat': stat,
        # 'stat_on': stat_on,
        'max_on': max_on,
        'cvi': paras_best['cvi'],
        'nsnv_use': paras_best['nsnv_use'],
    }

    statds = [statd | statd_para for statd in statds]

    return statds


def sex_of_phe(phe, phes_female):
    if phe in phes_female:
        sex = 'female'
    else:
        sex = 'both'
    return sex


# pheno=load from fphe
# phe=phenotype name
def load_best_or(
        dresult,
        methodms,
        method_modeld,
        model_integrate,
        fpheno,
        phes,
        phes_female,
        fcv,
        cvn,
        fymld,
        max_on,
        withcov_at_maxon,
        withcov,
        datasets,
        or_props):
    regond = load_regon(fymld)
    cv = sampleio.load_cv(fcv)
    pheno = sampleio.load_pheno(fpheno)

    stats = []
    for methodm in methodms:
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)

        for phe in phes:
            sex = sex_of_phe(phe, phes_female)

            paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                                va_on, fymld['cand'], method_modeld, None, model_integrate, cvn,)

            # print('paras_best_cvs', paras_best_cvs)

            for cvi in range(cvn):
                paras_best = paras_best_cvs[cvi]
                if paras_best is None:
                    continue

                method, mid_path = get_method_model_para(methodm, method_modeld, paras_best)
                for dataset in datasets:
                    score = load_sample_score_para(dresult, method, paras_best, cv, pheno, phe, dataset, withcov, reg_on, sex, mid_path)
                    statds = calc_best_stat_or(score, methodm, paras_best, dataset, phe, max_on, or_props)
                    stats += statds

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)
    return stats


def load_sample_score_para(dresult, method, paras_best, cv, pheno, phe, dataset, withcov, reg_on, sex, mid_path):
    fscore = iof.file_score(dresult, method, phe, paras_best['cvi'], paras_best, withcov, reg_on, mid_path)
    score = load_score(fscore, paras_best, cv, pheno, phe, dataset, sex)
    return score


def get_method_model_para(methodm, method_modeld, paras_best):
    # if methodm == 'boosting_integrate':
    if methodm.startswith('boosting_integrate'):
        # 'boosting_add' or 'boosting_nonadd'
        method_integrate = paras_best['method_integrate']
        mid_path = method_modeld[method_integrate]
    else:
        # None if not boosting
        mid_path = method_modeld[methodm]

    # 'boosting_nonadd' -> 'boosting'
    method = methodm_to_method(methodm)
    return (method, mid_path)


def load_sample_scored(
        dresult,
        methodms,
        method_modeld,
        model_integrate,
        fpheno,
        phe,
        sex,
        fcv,
        cvi,
        fymld,
        max_on,
        withcov_at_maxon,
        withcov,
        datasets):
    regond = load_regon(fymld)
    cv = sampleio.load_cv(fcv)
    pheno = sampleio.load_pheno(fpheno)

    # scored = {}
    # scored = defaultdict(lambda: {})
    # TODO: Is this all right?
    def recursive_defaultdict():
        return defaultdict(recursive_defaultdict)
    scored = recursive_defaultdict()
    for methodm in methodms:
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
        paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on, va_on,
                                            fymld['cand'], method_modeld, model_integrate, cvis=[cvi])
        paras_best = paras_best_cvs[cvi]
        if paras_best is None:
            continue

        method, mid_path = get_method_model_para(methodm, method_modeld, paras_best)

        for dataset in datasets:
            score = load_sample_score_para(dresult, method, paras_best, cv, pheno, phe, dataset, withcov, reg_on, sex, mid_path)
            scored[methodm][dataset] = score

    return scored


def extract_nsnv_best_para(paras_best_cvs, cvn, methodm, phe, max_on):
    stats = []
    for cvi in range(cvn):
        paras_best = paras_best_cvs[cvi]
        if paras_best is None:
            statd = None
        else:
            assert paras_best['cvi'] == cvi
            # print('paras_best', paras_best)
            statd = {
                # store 'boosting_nonadd' not 'boosting'
                'method': methodm,
                'phe': phe,
                'max_on': max_on,
                'cvi': paras_best['cvi'],
                # 'n': paras_best['n'],
                'nsnv_use': paras_best['nsnv_use'],
            }
        stats.append(statd)
    return stats


def load_best_nsnv(dresult, methodms, method_modeld, methodd_mode_va, model_integrate, phes, cvn, fymld, max_on, withcov_at_maxon):
    regond = load_regon(fymld)
    stats = []
    for methodm in methodms:
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
        for phe in phes:
            paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                                va_on, fymld['cand'], method_modeld, methodd_mode_va, model_integrate, cvn)
            nsnv_best_cvs = extract_nsnv_best_para(paras_best_cvs, cvn, methodm, phe, max_on)
            stats.extend(nsnv_best_cvs)

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)
    return stats


def load_best_para_table(dresult, methodms, method_modeld, model_integrate, phes, cvn, fymld, max_on, withcov_at_maxon):
    regond = load_regon(fymld)
    stats = []
    for methodm in methodms:
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
        for phe in phes:
            paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                                va_on, fymld['cand'], method_modeld, model_integrate, cvn)
            # nsnv_best_cvs = extract_nsnv_best_para(paras_best_cvs, cvn, methodm, phe, max_on)
            stats.extend(nsnv_best_cvs)

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)
    return stats


def load_para_grid(dresult, methodm, method_modeld, model_integrate, phe, cvi, fymld, max_on, plot_on,
                   withcov_at_maxon, withcov, dataset):
    regond = load_regon(fymld)

    reg_on, va_on = regon_vaon(methodm, regond, withcov)
    # reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)

    if methodm.startswith('boosting_integrate'):
        raise NotImplementedError('Should not use boosting_integrate.')

    method = methodm_to_method(methodm)
    mode = None
    # allstat for maxon
    allstat = load_allstat(dresult, method, phe, withcov, reg_on, method_modeld[methodm], mode)
    # logger.debug('allstat: {}'.format(allstat.head()))
    # logger.debug('allstat: {}'.format(allstat.info()))

    allstat_use = extract_cand_allstat(allstat, cvi, method, dataset, fymld['cand'],)

    # allstat_use = allstat.loc[(allstat['dataset'] == dataset) & (allstat['cvi'] == cvi), :]

    grid = allstat_use.pivot(index='nsnv_use', columns='lr', values=plot_on)

    # logger.debug('grid: {}'.format(grid))
    logger.debug('grid: {}'.format(grid.head()))

    # paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
    #                                    va_on, fymld['cand'], method_modeld, model_integrate, cvn)

    # nsnv_best_cvs = extract_nsnv_best_para(paras_best_cvs, cvn, methodm, phe, max_on)
    # stats.extend(nsnv_best_cvs)

    # remove None
    # stats = [x for x in stats if x is not None]
    # stats = pd.DataFrame.from_records(stats)

    parad_best = load_best_para(allstat, cvi, method, max_on, va_on, fymld['cand'])
    logger.debug('parad_best: {}'.format(parad_best))

    statd = {'grid': grid, 'allstat': allstat_use,
             'best_para': parad_best}
    return statd
    # return grid


def gmodels():
    models = ['overrec', 'rec', 'add', 'dom', 'overdom']
    return models


def gmodels_ss():
    gmodels = ['add', 'dom', 'rec', 'hetonly']
    return gmodels


def genetic_model_from_wgt(wgt):
    models = gmodels()

    # dfs = []
    # for phe_i in range(len(wgts)):
    uniq, counts = np.unique(wgt['boost_model_est'].dropna(), return_counts=True)
    d = dict(zip(uniq, counts))

    for k in models:
        if k not in d:
            d[k] = 0

    return d


def calc_genetic_model_para_best(dresult, method, methodm, paras_best, phe, cvi, max_on, mid_path):
    fwgt = iof.file_wgt(dresult, method, phe, cvi, paras_best, mid_path)
    wgt = wgtop.load_wgt_method(fwgt, method)
    # if fix_nsnv is None:
    #   nsnv=fix_nsnv
    # nsnv = paras_best['n']
    nsnv = paras_best['nsnv_use']

    wgt_nsnv = wgtop.format_wgt_boost(wgt, nsnv, '1')

    statd = genetic_model_from_wgt(wgt_nsnv)

    statd_para = {
        'method': methodm,
        'phe': phe,
        # 'max_on': max_on,
        'cvi': paras_best['cvi'],
        'nsnv_use': paras_best['nsnv_use'],
    }

    statd |= statd_para
    return statd


# def judge_genetic_model_gwas(ss):
#
#    # ss=ss.copy()
#
#    # col: '-log10(P)_add', dom, rec, (het)
#    # to distinguish overdom and overrec, effect_overdom is necessary
#    col_add = '-log10(P)_add'
#    col_dom = '-log10(P)_dom'
#    col_rec = '-log10(P)_rec'
#    col_overdom = '-log10(P)_overdom'
#    col_eff_overdom = 'effect_overdom'
#    col_model = 'model_est_gwas'
#
#    logger.debug('ss {}'.format(ss.head()))
#    ss.loc[:, col_model] = 'NaN'
#    if '-log10(P)_overdom' not in ss.columns:
#        raise RuntimeError('use het')
#        filt_add = (ss[col_add] > ss[col_dom]) & (ss[col_add] > ss[col_rec])
#        ss.loc[filt_add, col_model] = 'add'
#        filt_dom = (~filt_add) & (ss[col_dom] > ss[col_rec])
#        ss.loc[filt_dom, col_model] = 'dom'
#        filt_rec = (~filt_add) & (~filt_dom)
#        ss.loc[filt_rec, col_model] = 'rec'
#        assert (ss[col_model] != 'NaN').any()
#    else:
#        # do not forget to judge overdom or overrec
#        # assume p of overdom = overrec
#        # assume effect of overdom = overrec
#        filt_add = (ss[col_add] > ss[col_dom]) & (ss[col_add] > ss[col_rec]) & (ss[col_add] > ss[col_overdom])
#        assert (ss.loc[filt_add, col_model] == 'NaN').all()
#        ss.loc[filt_add, col_model] = 'add'
#
#        filt_dom = (~filt_add) & (ss[col_dom] > ss[col_rec]) & (ss[col_dom] > ss[col_overdom])
#        assert (ss.loc[filt_dom, col_model] == 'NaN').all()
#        ss.loc[filt_dom, col_model] = 'dom'
#
#        filt_rec = (~filt_add) & (~filt_dom) & (ss[col_rec] > ss[col_overdom])
#        assert (ss.loc[filt_rec, col_model] == 'NaN').all()
#        ss.loc[filt_rec, col_model] = 'rec'
#
#        # filt_overdom = (~filt_add) & (~filt_dom) & (~filt_rec) & (ss['effect_overdom'] > 0)
#        filt_overdom = (~filt_add) & (~filt_dom) & (~filt_rec) & (ss[col_eff_overdom] >= 0)
#        assert (ss.loc[filt_overdom, col_model] == 'NaN').all()
#        ss.loc[filt_overdom, col_model] = 'overdom'
#
#        filt_overrec = (~filt_add) & (~filt_dom) & (~filt_rec) & (ss[col_eff_overdom] < 0)
#        # filt_overrec = (~filt_add) & (~filt_dom) & (~filt_rec) & (~filt_overdom)
#        assert (ss.loc[filt_overrec, col_model] == 'NaN').all()
#        ss.loc[filt_overrec, col_model] = 'overrec'
#
#        # if one of -log10 is np.nan, this might happen
#        assert (ss[col_model] != 'NaN').any()
#
#    return ss


# def _merge_wgt_ss(wgt_ss, ss, gmodel, how='inner', use_sid_tmp=False):
#
#    if use_sid_tmp:
#        ss = ss.drop(columns=['id'])
#        ss = ss.rename(columns=lambda x: x if x == 'sid' else x + '_' + gmodel)
#        wgt_ss = pd.merge(wgt_ss, ss, on='sid', how=how)
#    else:
#        # avoid duplicated 'sid'
#        ss = ss.drop(columns=['sid'])
#        ss = ss.rename(columns=lambda x: x if x == 'id' else x + '_' + gmodel)
#        wgt_ss = pd.merge(wgt_ss, ss, on='id', how=how)
#
#    # logger.debug('ss {}'.format(ss.head()))
#    # logger.debug('ss {}'.format(ss.columns))
#
#    # avoid contamination of ss.columns
#    # ss = ss.copy()
#    # ss.columns = [x + '_' + gmodel for x in ss.columns]
#    # logger.debug('ss {}'.format(ss.head()))
#    # if use_sid_tmp:
#    #    wgt_ss = pd.merge(wgt_ss, ss, left_on='sid', right_on='sid_' + gmodel, how=how)
#    #    # wgt_ss = pd.merge(wgt_ss, ss, left_on='sid', right_on='sid_' + gmodel, how='inner')
#    # else:
#    #    wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how=how)
#    #    # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how='inner')
#
#    # wgt_ss = wgt_ss.drop(columns=['id_' + gmodel, 'sid_' + gmodel])
#    # wgt_ss=wgt_ss.drop(columns=['id_'+gmodel])
#
#    return wgt_ss


# def merge_wgt_ssd(wgt_ss, ssd, how='inner', use_sid_tmp=False):
#    # Does not assume ss-only variants exist
#    # TMP
#    ss_cols_use = ['id', 'sid', '-log10(P)', 'effect', 'p']
#    # ss_cols_use = ['id', '-log10(P)', 'effect', 'p']
#    for gmodel, ss in ssd.items():
#        ss = ss.loc[:, ss_cols_use].copy()
#        logger.debug('ss {}'.format(ss.head()))
#
#        if gmodel == 'hetonly':
#            for gmodel_merge in ['overdom', 'overrec']:
#                wgt_ss = _merge_wgt_ss(wgt_ss, ss, gmodel_merge, how, use_sid_tmp)
#                # wgt_ss = merge_wgt_ss(wgt_ss, ss, gmodel_merge, use_sid_tmp)
#
#            # ss.columns = [x + '_overdom' for x in ss_cols_use]
#            # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_overdom', how='inner')
#            # ss.columns = [x + '_overrec' for x in ss_cols_use]
#            # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_overrec', how='inner')
#        else:
#            wgt_ss = _merge_wgt_ss(wgt_ss, ss, gmodel, how, use_sid_tmp)
#
#            # [ss.columns = [x + '_' + gmodel for x in ss_cols_use]
#            # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how='inner')
#        logger.debug('merge wgt_ss {}'.format(wgt_ss.head()))
#    return wgt_ss


def calc_genetic_model_gwas_para_best(dresult, method, methodm, paras_best, phe, cvi, max_on, mid_path):

    use_sid_tmp = ('dataset.12' in str(dresult))

    fwgt = iof.file_wgt(dresult, method, phe, cvi, paras_best, mid_path)
    wgt = wgtop.load_wgt_method(fwgt, method)
    # if fix_nsnv is None:
    #   nsnv=fix_nsnv
    # nsnv = paras_best['n']
    nsnv = paras_best['nsnv_use']

    wgt_nsnv = wgtop.format_wgt_boost(wgt, nsnv, '1')

    # load ss
    # for gmodel in gmodels:
    #    fss=load_ss()

    gmodels = gmodels_ss()
    ssd = ssop.load_ssd_raw(dresult, phe, cvi, gmodels)

    wgt_ss = wgt_nsnv
    logger.debug('wgt_ss {}'.format(wgt_ss.head()))

    wgt_ss = ssop.merge_wgt_ssd(wgt_ss, ssd, use_sid_tmp=use_sid_tmp)

    # TMP
    # ss_cols_use = ['id', 'sid', '-log10(P)', 'effect', 'p']
    # ss_cols_use = ['id', '-log10(P)', 'effect', 'p']
    # for gmodel, ss in ssd.items():
    #    ss = ss.loc[:, ss_cols_use].copy()
    #    logger.debug('ss {}'.format(ss.head()))

    #    if gmodel == 'hetonly':
    #        for gmodel_merge in ['overdom', 'overrec']:
    #            wgt_ss = merge_wgt_ss(wgt_ss, ss, gmodel_merge, use_sid_tmp)
    #        # ss.columns = [x + '_overdom' for x in ss_cols_use]
    #        # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_overdom', how='inner')
    #        # ss.columns = [x + '_overrec' for x in ss_cols_use]
    #        # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_overrec', how='inner')
    #    else:
    #        wgt_ss = merge_wgt_ss(wgt_ss, ss, gmodel, use_sid_tmp)
    #        # [ss.columns = [x + '_' + gmodel for x in ss_cols_use]
    #        # wgt_ss = pd.merge(wgt_ss, ss, left_on='id', right_on='id_' + gmodel, how='inner')
    #    logger.debug('merge wgt_ss {}'.format(wgt_ss.head()))

    wgt_ss = ssop.judge_genetic_model_gwas(wgt_ss)

    # statd = genetic_model_from_wgt(wgt_nsnv)

    # statd = {
    #    'method': methodm,
    #    'phe': phe,
    #    'max_on': max_on,
    #    'cvi': paras_best['cvi'],
    #    'nsnv_use': paras_best['nsnv_use'],
    # }

    # statd['wgt'] = wgt
    return wgt_ss


# def load_best_genetic_model(dresult, methodm, paras_best_cvs, phe, max_on, cvi, method_modeld):
#    # if methodm == 'boosting_integrate':
#    if methodm.startswith('boosting_integrate'):
#        raise RuntimeError('Do not use boosting_integrate here.')
#
#    # TMP
#    if methodm != 'boosting_nonadd':
#        raise RuntimeError('Use boosting_nonadd only.')
#
#    # 'boosting_nonadd' -> 'boosting'
#    method = methodm_to_method(methodm)
#
#    # None if not boosting
#    model = method_modeld[methodm]
#
#    stats = []
#    # for cvi in range(cvn):
#    paras_best = paras_best_cvs[cvi]
#    statd = calc_genetic_model_para_best(dresult, method, methodm, paras_best, phe, cvi, max_on, model)
#    stats.append(statd)
#    return stats


def load_genetic_model(dresult, method_modeld, phes, cvi, fymld, max_on, withcov_at_maxon):
    regond = load_regon(fymld)
    stats = []

    methodm = 'boosting_nonadd'
    reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
    for phe in phes:
        paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                            va_on, fymld['cand'], method_modeld, None, None, cvis=[cvi])

        # stats_best_cvs = load_best_genetic_model(dresult, methodm, paras_best_cvs, phe, max_on, cvi, method_modeld)

        # 'boosting_nonadd' -> 'boosting'
        method = methodm_to_method(methodm)
        mid_path = method_modeld[methodm]

        paras_best = paras_best_cvs[cvi]
        statd = calc_genetic_model_para_best(dresult, method, methodm, paras_best, phe, cvi, max_on, mid_path)
        stats.append(statd)

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)
    return stats


def load_genetic_model_gwas(dresult, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon):
    regond = load_regon(fymld)
    # stats = []
    # statsd = {}

    methodm = 'boosting_nonadd'
    reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)

    # for phe in phes:
    paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                        va_on, fymld['cand'], method_modeld, None,None, cvis=[cvi])

    # stats_best_cvs = load_best_genetic_model(dresult, methodm, paras_best_cvs, phe, max_on, cvi, method_modeld)

    # 'boosting_nonadd' -> 'boosting'
    method = methodm_to_method(methodm)
    mid_path = method_modeld[methodm]

    paras_best = paras_best_cvs[cvi]
    stat = calc_genetic_model_gwas_para_best(dresult, method, methodm, paras_best, phe, cvi, max_on, mid_path)
    # TODO: to_csv(statd['wgt'])
    # stats.append(statd)
    # statsd[phe] = statd

    # remove None
    # stats = [x for x in stats if x is not None]
    # stats = pd.DataFrame.from_records(stats)
    # return statsd
    return stat


# def load_genetic_model_vs_nsnv(dresult, method_modeld, phes, cvi, fymld, max_on, withcov_at_maxon):
#    regond = load_regon(fymld)
#    stats = []
#
#    methodm = 'boosting_nonadd'
#    reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
#    for phe in phes:
#        paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
#                                            va_on, fymld['cand'], method_modeld, None,None, cvis=[cvi])
#
#        # stats_best_cvs = load_best_genetic_model(dresult, methodm, paras_best_cvs, phe, max_on, cvi, method_modeld)
#
#        # 'boosting_nonadd' -> 'boosting'
#        method = methodm_to_method(methodm)
#        mid_path = method_modeld[methodm]
#
#        paras_best = paras_best_cvs[cvi]
#        statd = calc_genetic_model_para_best(dresult, method, methodm, paras_best, phe, cvi, max_on, mid_path)
#
#        # FIXME:  add nsnvs
#        statd['nsnv_use'] = paras_best['nsnv_use']
#
#        stats.append(statd)
#
#    # remove None
#    stats = [x for x in stats if x is not None]
#    stats = pd.DataFrame.from_records(stats)
#    return stats


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
# def align_a1_to_alt_ss(ss):
#    filt_unmatch = (ss['alt'] != ss['a1'])
#
#    ss.loc[filt_unmatch, 'a1_tmp'] = ss.loc[filt_unmatch, 'a2']
#    ss.loc[filt_unmatch, 'a2'] = ss.loc[filt_unmatch, 'a1']
#    ss.loc[filt_unmatch, 'a1'] = ss.loc[filt_unmatch, 'a1_tmp']
#
#    assert (ss['alt'] == ss['a1']).all()
#    assert (ss['ref'] == ss['a2']).all()
#
#    ss.loc[filt_unmatch, 'effect'] = -ss.loc[filt_unmatch, 'effect']
#    if 'or' in ss.columns:
#        ss.loc[filt_unmatch, 'or'] = 1.0 / ss.loc[filt_unmatch, 'or']
#    ss.loc[filt_unmatch, 'l95'] = 1.0 / ss.loc[filt_unmatch, 'u95']
#
#    return ss


# def load_ss(fss, align_a1_alt=False, maf_thre=None):
#     '''
#     id: rsid
#     sid: 'chrom':'pos'
#     sida: 'chrom':'pos':'a1':'a2'
#     '''
#     # TMP; for condrec, condadd of cv1-
#     if not os.path.exists(fss):
#         logger.info("skip since fss does not exist: {}".format(fss))
#         return None

#     ss = pd.read_csv(fss, sep='\t', dtype={'P': str})
#     ss = ss.rename(columns=str.lower)

#     if align_a1_alt:
#         ss = align_a1_to_alt_ss(ss)

#     ss = snvop.add_sids(ss)

#     ss.loc[:, 'effect'] = np.log(ss['or'])

#     # ss = ss.rename(columns={'a1': 'A1', 'a2': 'A2'})
#     # TMP
#     # if 'dataset.12' in str(fss):
#     # ss.loc[:, 'id'] = ss['sida']

#     # def str_from_log10_decimal(x):
#     #    # https://stackoverflow.com/questions/36092203/why-is-1e400-not-an-int
#     #    dec = Decimal(10)**(Decimal(-1) * Decimal(x))
#     #    if dec < 10**(-4):
#     #        num_str = '{:.6e}'.format(dec)
#     #    else:
#     #        num_str = '{:.6f}'.format(dec)
#     #    return num_str

#     def p_decimal_from_str(x):
#         # https://stackoverflow.com/questions/36092203/why-is-1e400-not-an-int
#         return Decimal(x)

#     def log10_p_from_str(x):
#         return float(Decimal(x).log10())

#     # error
#     # def p_f128_from_str(x):
#     #    return np.float128(Decimal(x))

#     ss.loc[:, 'p_str'] = ss['p']

#     # float
#     # np.float128 prints out 0.0 but actually 1e-334. confirmed by np.log10
#     ss.loc[:, 'p'] = ss['p_str'].astype(np.float128)
#     ss.loc[:, '-log10(P)'] = -ss['p'].map(np.log10).astype(np.float64)

#     if (ss['p'] == 0.0).any():
#         logger.error('log10(P): {}'.format((ss['p'] == 0.0).any()))
#         logger.error('log10(P): {}'.format(ss.loc[ss['p'] == 0.0, 'p']))
#         raise RuntimeError

#     # check if p and log10(p) work fine
#     # if -log10(p) and decimal match, it's fine
#     ss_check = ss[['p_str', 'p', '-log10(P)']].copy()
#     ss_check.loc[:, 'p_f64'] = ss_check['p_str'].astype(np.float64)
#     # decimal from str: actually you can run log from decimal 'p'
#     ss_check.loc[:, '-log10(P)_decimal'] = -ss_check['p_str'].apply(lambda x: log10_p_from_str(x))
#     # decimal type
#     ss_check.loc[:, 'p_decimal'] = ss_check['p_str'].apply(lambda x: p_decimal_from_str(x))
#     assert np.allclose(ss_check['-log10(P)'], ss_check['-log10(P)_decimal'])
#     if (ss_check['p_f64'] == 0.0).any():
#         logger.info('p_f64: {}'.format(ss_check.loc[ss_check['p_f64'] == 0.0, ['p', 'p_f64', 'p_decimal', '-log10(P)', '-log10(P)_decimal']]))

#     if maf_thre is not None:
#         raise NotImplementedError
#         # raise RuntimeError('ny for dataset<=11; not necessary for dataset12')
#         # TODO: use ukb_imp_com.frq for maf
#         # print('limit maf', maf_thre)
#         # fco_phe = fco.replace('%', phenotype)
#         # gt = pd.read_hdf(fco_phe + '.obs.gt', key='gt')
#         # ss = pd.merge(ss, gt[['vid', 'maf']], left_on='sida', right_on='vid', how='left')
#         # ss = ss[ss['maf'] >= maf_thre].copy()

#     logger.debug('ss len: {}'.format(len(ss)))

#     return ss


def extract_ss_paper(ss):
    # ref,alt
    cols = 'chrom,pos,id,a1,a2,or,p,-log10(P),effect,sid,sida'.split(',')
    ss = ss.loc[:, cols].copy()
    return ss


def search_iter_loss(dresult, phe, cvi, lr, mid_path):
    # find the minimim iter of existing file
    for iteri in range(30):
        floss = iof.file_loss(dresult, phe, cvi, lr, iteri, mid_path)
        print('floss', floss)
        if os.path.isfile(floss):
            return iteri
    raise RuntimeError('loss does not exist: ', phe, cvi, lr, mid_path)


def load_loss(floss):
    loss = pd.read_csv(floss, sep='\t', header=0, dtype={'loss': np.float64})
    loss = loss.rename(columns=str.lower)

    logger.debug('loss {}'.format(loss))

    loss = snvop.add_sids(loss)
    loss.loc[:, 'id'] = loss['rs']

    # due to rust loss
    if (loss['loss'] > 1.0 * 10**300).any():
        loss_max = loss.loc[loss['loss'] < 1.0 * 10**300, 'loss'].max()
        logger.debug('loss_max: {}'.format(loss_max))
        loss.loc[loss['loss'] > 1.0 * 10**300, 'loss'] = loss_max
    return loss


# def load_ssd_raw(dresult, phe, cvi, gmodels=['add', 'dom', 'rec', 'hetonly'], mid_path=None):
#    # gwas
#    ssd = {}
#    # gmodels = ['add', 'dom', 'rec', 'hetonly']
#    for gmodel in gmodels:
#        fss = iof.file_ss(dresult, phe, cvi, gmodel, mid_path)
#        ss = ssio.load_ss(fss)
#        # ss = extract_ss_paper(ss)
#        ssd[gmodel] = ss
#    return ssd


def update_ssd_wgt(ssd, wgt, methodm, use_sid_tmp=False):
    # 'boosting_nonadd' -> 'boosting'
    # method = methodm_to_method(methodm)

    if wgt is None:
        logger.warning('wgt is None: {}.'.format(methodm))
        return ssd

    for gmodel in ssd:
        ss = ssd[gmodel]
        if methodm == 'boosting_loss':
            if use_sid_tmp:
                ss = pd.merge(ss, wgt[['sid', 'loss']], on='sid', how='left')
            else:
                ss = pd.merge(ss, wgt[['id', 'loss']], on='id', how='left')

            # FIXME: no need to align a1??
            # -> should check already aligned as in 'boosting'
            # -> should be no problem for plotting

            logger.debug('loss max: {}'.format(ss['loss'].max()))
            ss.loc[:, 'decreasing_loss'] = -ss['loss']

        # elif method == 'boosting':
        elif methodm == 'boosting_nonadd':
            # if you want to update col, also update cols_use below
            cols = 'id,sid,a1,a2,boost_score1_rel,boost_score2_rel,boost_model_est'.split(',')
            # cols = 'id,sid,a1,a2,boost_score1_rel,boost_score2_rel'.split(',')
            wgt_tmp = wgt[cols].copy()
            # wgt_tmp = wgt.copy()

            # align a1, a2
            if use_sid_tmp:
                wgt_tmp = pd.merge(wgt_tmp, ss[['sid', 'a1', 'a2']], on='sid', how='left', suffixes=('', '_ss'))
            else:
                wgt_tmp = pd.merge(wgt_tmp, ss[['id', 'a1', 'a2']], on='id', how='left', suffixes=('', '_ss'))
            if (wgt_tmp['a1'] != wgt_tmp['a1_ss']).any() and (wgt_tmp['a2'] != wgt_tmp['a2_ss']).any():
                raise NotImplementedError('ny: troublesome for boosting.')

            # logger.debug('wgt_tmp here {}'.format(wgt_tmp))

            wgt_tmp.loc[:, methodm + '_effect'] = wgt_tmp['boost_score1_rel']
            # wgt[method + '_effect'] = wgt['boost_score1'] - wgt['boost_score0']
            wgt_tmp.loc[:, methodm + '_effect_abs'] = np.abs(wgt_tmp[methodm + '_effect'])

            wgt_tmp.loc[:, methodm + '_effect2'] = wgt_tmp['boost_score2_rel']
            wgt_tmp.loc[:, methodm + '_effect2_abs'] = np.abs(wgt_tmp[methodm + '_effect2'])

            cols_use = 'boost_model_est'.split(',')
            cols_use_effect = 'effect,effect_abs,effect2,effect2_abs'.split(',')
            cols_use += [methodm + '_' + x for x in cols_use_effect]

            # logger.debug('wgt_tmp here {}'.format(wgt_tmp))

            # if 'dataset.12' in str(fss):
            if use_sid_tmp:
                # TMP
                ss = pd.merge(ss, wgt_tmp[['sid'] + cols_use], on='sid', how='left')
                # ss = pd.merge(ss, wgt_tmp[['sid', methodm + '_effect']], on='sid', how='left')
                # ss = pd.merge(ss, wgt_tmp[[method + '_sid', method + '_effect']], left_on='sid', right_on=method + '_sid', how='left')
            else:
                ss = pd.merge(ss, wgt_tmp[['id'] + cols_use], on='id', how='left')
                # ss = pd.merge(ss, wgt_tmp[['id', methodm + '_effect']], on='id', how='left')

        else:
            # additive model
            # 'boosting_add' or prev methods

            logger.debug('methodm {}'.format(methodm))
            logger.debug('wgt {}'.format(wgt.columns))
            cols = 'id,sid,a1,a2,effect'.split(',')
            wgt_tmp = wgt[cols].copy()

            # align a1, a2
            if use_sid_tmp:
                wgt_tmp = pd.merge(wgt_tmp, ss[['sid', 'a1', 'a2']], on='sid', how='left', suffixes=('', '_ss'))
            else:
                wgt_tmp = pd.merge(wgt_tmp, ss[['id', 'a1', 'a2']], on='id', how='left', suffixes=('', '_ss'))
            filt_flip = wgt_tmp['a1'] != wgt_tmp['a1_ss']
            assert (wgt_tmp.loc[filt_flip, 'a1'] == wgt_tmp.loc[filt_flip, 'a2_ss']).all()
            wgt_tmp.loc[filt_flip, 'effect'] = - wgt_tmp.loc[filt_flip, 'effect']

            wgt_tmp.loc[:, methodm + '_effect'] = wgt_tmp['effect']
            wgt_tmp.loc[:, methodm + '_effect_abs'] = np.abs(wgt_tmp[methodm + '_effect'])

            cols_use = 'effect,effect_abs'.split(',')
            cols_use = [methodm + '_' + x for x in cols_use]

            if use_sid_tmp:
                ss = pd.merge(ss, wgt_tmp[['sid'] + cols_use], on='sid', how='left')
            else:
                ss = pd.merge(ss, wgt_tmp[['id'] + cols_use], on='id', how='left')

        logger.debug('ss {}'.format(ss))
        ssd[gmodel] = ss
    return ssd


def load_best_wgt(dresult, methodm, phe, cvi, paras_best_cvs, method_modeld):
    # if methodm == 'boosting_integrate':
    if methodm.startswith('boosting_integrate'):
        # TODO: implement code are already in public_wgt.py:load_best_wgt()
        raise RuntimeError('Do not use boosting_integrate here.')

    # 'boosting_nonadd' -> 'boosting'
    method = methodm_to_method(methodm)
    # None if not boosting
    mid_path = method_modeld[methodm]

    if methodm == 'boosting_loss':
        # loss
        lr_loss = '1.0'
        iter_loss = search_iter_loss(dresult, phe, cvi, lr_loss, mid_path)
        floss = iof.file_loss(dresult, phe, cvi, lr_loss, iter_loss, mid_path)
        loss = load_loss(floss)
        wgt = loss
    elif method == 'boosting':
        # elif methodm == 'boosting_nonadd':
        # if (methodm != 'boosting_nonadd'):
        # raise RuntimeError('Use boosting_nonadd only: {}'.format(methodm))
        paras_best = paras_best_cvs[cvi]
        logger.debug('paras_best method phe: {}, {}, {}'.format(methodm, phe, paras_best))
        # should be methodm='boosting_nonadd'
        fwgt = iof.file_wgt(dresult, method, phe, cvi, paras_best, mid_path)
        wgt = wgtop.load_wgt_method(fwgt, method)
        # if fix_nsnv is None:
        #   nsnv=fix_nsnv
        # nsnv = paras_best['n']
        nsnv = paras_best['nsnv_use']
        logger.debug('nsnv {}'.format(nsnv))

        wgt = wgtop.format_wgt_boost(wgt, nsnv, '1')

        # logger.debug('wgt {}',wgt.columns)
    else:
        # boosting_add or prev methods
        paras_best = paras_best_cvs[cvi]
        fwgt = iof.file_wgt(dresult, method, phe, cvi, paras_best, mid_path)
        try:
            wgt = wgtop.load_wgt_method(fwgt, method)
        except OSError as e:
            logger.warning(e)
            wgt = None
    return wgt


def load_manhattan_ssd(dresult, methodms, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon, gmodels):
    '''
    Load ss containing loss, gwas of selected snvs in wgt
    '''
    regond = load_regon(fymld)
    # stats = []

    # ssd = load_ssd_raw(dresult, phe, cvi, ['add'])
    ssd = ssop.load_ssd_raw(dresult, phe, cvi, gmodels)

    methodms_use = []
    for methodm in methodms:
        # methodm = 'boosting_nonadd'
        reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
        paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                            va_on, fymld['cand'], method_modeld, None, None, cvis=[cvi])

        # stats_best_cvs = load_best_genetic_model(dresult, methodm, paras_best_cvs, phe, reg_on, max_on, cvn, method_modeld)
        # ssd = load_best_ssd(dresult, methodm, paras_best_cvs, phe, max_on, method_modeld, cvi)

        wgt = load_best_wgt(dresult, methodm, phe, cvi, paras_best_cvs, method_modeld)
        if wgt is None:
            continue
        else:
            methodms_use.append(methodm)
        ssd = update_ssd_wgt(ssd, wgt, methodm, use_sid_tmp=('dataset.12' in str(dresult)))

    # for loss
    methodm = 'boosting_loss'
    # methodms_use.append(methodm)
    wgt = load_best_wgt(dresult, methodm, phe, cvi, None, method_modeld)
    ssd = update_ssd_wgt(ssd, wgt, methodm, use_sid_tmp=('dataset.12' in str(dresult)))

    # stats.extend(stats_best_cvs)

    # remove None
    # stats = [x for x in stats if x is not None]
    # stats = pd.DataFrame.from_records(stats)
    return ssd, methodms_use


# plot_on for auc
def get_withcovs(max_on, plot_on=None,
                 withcov_only=False, nocov_only=False,
                 ):
    # since sometimes needs
    # [(F,F), (F,T), (T,T)]
    # especially for auc
    #
    #
    # if maxon=nagel and ploton=nagel
    # [(T,T)]
    #
    # if maxon=nagel
    # [(F,T),(T,T)]
    #
    # if ploton=nagel
    # [(T,F),(T,T)]
    #
    # else
    # [(F,F),(F,T),(T,F),(T,T)]

    # (withcov, withcov_at_maxon)

    if max_on == 'nagelkerke' and plot_on == 'nagelkerke':
        # nocovs_and_maxon=[(True, True),(True, False)]
        withcovs = [(True, True)]
        # return [(True, True)]

    elif max_on == 'nagelkerke':
        withcovs = [(False, True), (True, True)]
        # return [(False, True), (True, True)]

    elif plot_on == 'nagelkerke':
        withcovs = [(True, False), (True, True)]

    else:
        withcovs = [(False, False), (False, True), (True, False), (True, True)]

    if withcov_only and nocov_only:
        raise RuntimeError('Wrong arg')
    elif withcov_only:
        return [x for x in withcovs if x[0]]
    elif nocov_only:
        return [x for x in withcovs if not x[0]]
    else:
        return withcovs


def get_withcov_at_maxon(max_on):
    # withcov_at_maxon only
    # for nsnv

    if max_on == 'nagelkerke':
        return [True]

    return [False, True]


def dname_maxon(droot, max_on=None, withcov=None, withcov_at_maxon=None, plot_on=None, phe=None, cvi=None, dataset=None, method=None, sex=None):
    dl = []
    if max_on is not None:
        dl += ['maxon_' + max_on]
    if plot_on is not None:
        dl += ['ploton_' + plot_on]

    if withcov_at_maxon is None:
        pass
    elif withcov_at_maxon:
        dl += ['maxon_withcov']
    else:
        dl += ['maxon_nocov']

    if withcov is None:
        pass
    elif withcov:
        dl += ['withcov']
    else:
        dl += ['nocov']
    d = '.'.join(dl)

    # d = 'maxon_' + max_on
    # if withcov_at_maxon is None:
    #    pass
    # elif withcov_at_maxon:
    #    d += '.maxon_withcov'
    # else:
    #    d += '.maxon_nocov'

    # if withcov is None:
    #    pass
    # elif withcov:
    #    d += '.withcov'
    # else:
    #    d += '.nocov'

    d2l = []
    if method is not None:
        d2l += [method]
    if phe is not None:
        d2l += [phe]
    if cvi is not None:
        d2l += ['cv' + str(cvi)]
    if dataset is not None:
        d2l += [dataset]
    if sex is not None:
        d2l += [sex]
    d2 = '.'.join(d2l)

    # if d='', d=droot or droot/d2
    # if d2='', d=droot/d
    d = droot / d / d2

    # if phe is None and cvi is None:
    #    d = droot / d
    # else:
    #    if phe is not None and cvi is None:
    #       d2 = phe
    #    elif phe is None and cvi is not None:
    #       d2 = 'cv' + str(cvi)
    #    else:
    #       d2 = phe + '.cv' + str(cvi)
    #    d = droot / d / d2

    # good idea
    # ds=['maxon_'+max_on,'maxonnocov_'+str(nocov_at_maxon),'nocov_'+str(nocov)]
    # dname='.'.join(ds)
    # d=droot+dname+'/'
    return d


# def dname_maxon_create(droot, max_on, withcov=None, withcov_at_maxon=None, phe=None, cvi=None,dataset=None):
def dname_maxon_create(*args, **kwargs):
    d = dname_maxon(*args, **kwargs)
    # d = dname_maxon(droot, max_on, withcov, withcov_at_maxon, phe, cvi)
    os.makedirs(d, exist_ok=True)
    return d


def run_best_para_table(droot, pformats, dresult, methodms, method_modeld, model_integrate, phes, cvn, fymld):
    datasets = ['ts']

    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))

            stats = load_best_para_table(dresult, methodms, method_modeld, model_integrate, phes, cvn, fymld, max_on, withcov_at_maxon)
            logger.debug('stats_nsnv\n{}'.format(stats.head()))

            # dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon)
            # output to tsv

            # plot
            # dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon)
            # paper_validate.validate_nsnv(dval, pformats, stats, methodms, phes)


def run_para_grid(droot, pformats, dresult, methodms, method_modeld, model_integrate, phes, cvis, fymld):
    datasets = ['ts']

    for max_on in ['nagelkerke']:
        for plot_on in ['nagelkerke']:
            # for withcov in get_withcov_at_maxon(plot_on):
            withcovs = [(True, True)]
            for (withcov, withcov_at_maxon) in withcovs:
                logger.debug('withcov: {}'.format(withcov))

                for phe in phes:
                    for cvi in cvis:
                        for dataset in datasets:
                            for methodm in methodms:
                                stats = load_para_grid(dresult, methodm, method_modeld, model_integrate, phe,
                                                       cvi, fymld, max_on, plot_on, withcov_at_maxon, withcov, dataset)
                                logger.debug('stats grid\n{}'.format(stats['grid'].head()))

                                dval = dname_maxon_create(droot, max_on=max_on, withcov=withcov, plot_on=plot_on,
                                                          phe=phe, cvi=cvi, dataset=dataset, method=methodm)

                                # output to tsv
                                stats['grid'].to_csv(dval / 'grid.tsv', na_rep='NaN', sep='\t')

                                # plot
                                paper_validate.validate_para_grid(dval, pformats, stats, plot_on, methodm, phe)


# def convert_relative_stat(stats: pd.DataFrame, stats_covonly: pd.DataFrame):
#    baseline_method = 'boosting_integrate'
#    # cols_index = ['method', 'i_cv', 'nsnvs', 'phe']
#    cols_index = ['method', 'cvi', 'phe']
#    stats_out = stats.copy().set_index(cols_index)
#    stats_covonly = stats_covonly.copy().set_index(cols_index)
#    print('stats_out', stats_out.head(50))
#    print('stats_covonly', stats_covonly.head(50))
#    for phe in stats['phe'].unique():
#        for i_cv in stats['cvi'].unique():
#            auc_covonly = stats_covonly.loc[('covonly', i_cv, phe), 'stat']
#            base_auc = stats_out.loc[(baseline_method, i_cv, phe), 'stat']
#            # checked before
#            # assert base_auc >= auc_covonly
#            # skip if base_auc >= auc_covonly
#            if base_auc < auc_covonly:
#                print('baseline_method',  baseline_method)
#                print(' phe, icv', phe, i_cv)
#                print('covonly, base', auc_covonly, base_auc)
#                print("stat of baseline is smaller than covonly")
#                print("SKIP this icv")
#                continue
#            # this results in negative sign
#            # if base_auc < auc_covonly:
#            #    print('plot_on,baseline_method',plot_on,baseline_method)
#            #    print('nsnv, phe, icv',nsnv,phe,i_cv)
#            #    print('covonly, base',auc_covonly, base_auc)
#            #    raise RuntimeError("stat of baseline is smaller than covonly")
#            for method in stats['method'].unique():
#                index = (method, i_cv, phe)
#                if index not in stats_out.index:
#                    continue
#                auc_method = stats_out.loc[index, 'stat']
#                auc_new = (auc_method - base_auc) / (base_auc - auc_covonly)
#                # auc_new=100*(auc_method-base_auc)/(base_auc-auc_covonly)
#                stats_out.loc[index, 'stat'] = auc_new
#                # stats_out.loc[index,plot_on]=auc_method-base_auc


def convert_median_relative_stat(stats_median: pd.DataFrame, baseline_method='boosting_integrate',
                                 kind='relative'):
    # for nagelkerke only

    # baseline_method = 'boosting_integrate'
    # cols_index = ['method', 'i_cv', 'nsnvs', 'phe']
    cols_index = ['method', 'phe']
    phes = stats_median['phe'].unique()
    methods = stats_median['method'].unique()
    stats_median = stats_median.copy().set_index(cols_index)
    stats_out = stats_median.copy()
    #: clean df
    stats_out.loc[:, :] = np.nan
    # stats_covonly = stats_covonly.copy().set_index(cols_index)
    # print('stats_out', stats_out.head(50))
    # print('stats_covonly', stats_covonly.head(50))
    logger.debug('stats_out {}'.format(stats_out))
    for phe in phes:
        # for phe in stats_median['phe'].unique():
        # for i_cv in stats['cvi'].unique():
        # auc_covonly = stats_covonly.loc[('covonly', i_cv, phe), 'stat']

        if (baseline_method, phe) not in stats_out.index:
            logger.debug('skip baseline, phe: {}, {}'.format(baseline_method, phe))
            continue
        base_auc = stats_median.loc[(baseline_method, phe), 'stat']
        # base_auc = stats_out.loc[(baseline_method, phe), 'stat']

        # checked before
        # assert base_auc >= auc_covonly
        # skip if base_auc >= auc_covonly
        # if base_auc < auc_covonly:
        #    print('baseline_method', baseline_method)
        #    print(' phe, icv', phe, i_cv)
        #    print('covonly, base', auc_covonly, base_auc)
        #    print("stat of baseline is smaller than covonly")
        #    print("SKIP this icv")
        #    continue
        # this results in negative sign
        # if base_auc < auc_covonly:
        #    print('plot_on,baseline_method',plot_on,baseline_method)
        #    print('nsnv, phe, icv',nsnv,phe,i_cv)
        #    print('covonly, base',auc_covonly, base_auc)
        #    raise RuntimeError("stat of baseline is smaller than covonly")

        # for method in stats_median['method'].unique():
        for method in methods:
            index = (method, phe)
            if index not in stats_median.index:
                continue
            auc_method = stats_median.loc[index, 'stat']
            if kind == 'relative':
                auc_new = (auc_method - base_auc) / (base_auc)
            elif kind == 'diff':
                auc_new = (auc_method - base_auc)
            else:
                raise NotImplementedError
            # auc_new = (auc_method - base_auc) / (base_auc - auc_covonly)
            # auc_new=100*(auc_method-base_auc)/(base_auc-auc_covonly)
            stats_out.loc[index, 'stat'] = auc_new
            # stats_out.loc[index,plot_on]=auc_method-base_auc

    return stats_out


# def calculate_top_number(stats):
#    # cols_index=['method','nsnvs','phe']
#    # countd = defaultdict(lambda: 0)
#    countd = defaultdict(lambda: 0)
#    phed = defaultdict(lambda: [])
#    for phe in stats['phe'].unique():
#        logger.debug('phe: {}'.format(stats.loc[stats['phe'] == phe, 'stat'].idxmax()))
#        idx = stats.loc[stats['phe'] == phe, 'stat'].idxmax()
#        method_max = stats.loc[idx, 'method']
#        countd[method_max] += 1
#        phed[method_max] += [phe]
#    stats_top = pd.DataFrame([countd, phed], index=['top_count', 'phes'])
#    stats_top = stats_top.T
#    stats_top = stats_top.reset_index(drop=False)
#
#    stats_top.loc[:,'phes']=stats_top['phes'].apply(lambda x: ','.join(x))
#    #stats_rank.loc[:,'phes']=stats_rank['phes'].apply(lambda x: ','.join(x))
#
#    # stats_top = pd.DataFrame([countd,phed], columns=['top_count','phes'])
#    # stats_top = pd.DataFrame.from_dict(countd, orient='index', columns=['top_count'])
#    return stats_top


def calculate_rank(stats):
    # cols_index=['method','nsnvs','phe']
    # countd = defaultdict(lambda: 0)
    method_rank_phed = defaultdict(lambda: [])

    # countd = defaultdict(lambda: 0)
    # phed = defaultdict(lambda: [])
    for phe in stats['phe'].unique():
        logger.debug('phe: {}'.format(stats.loc[stats['phe'] == phe, 'stat'].idxmax()))
        stats_phe = stats.loc[stats['phe'] == phe, :].copy()
        stats_phe = stats_phe.sort_values('stat', ascending=False).reset_index(drop=True)
        for (method, rank) in zip(stats_phe['method'], stats_phe.index):
            # +1 since rank starts from 0
            method_rank_phed[(method, rank + 1)] += [phe]

    rankd = []
    for (method, rank) in method_rank_phed:
        rankd += [{'method': method, 'rank': rank, 'phes': method_rank_phed[(method, rank)]}]
    # print(method_rank_phed)
    # stats_rank = pd.DataFrame.from_dict(method_rank_phed)
    stats_rank = pd.DataFrame.from_dict(rankd)

    stats_rank.loc[:, 'num_phes'] = stats_rank['phes'].apply(lambda x: len(x))
    stats_rank.loc[:, 'phes'] = stats_rank['phes'].apply(lambda x: ', '.join(x))
    stats_rank = stats_rank.loc[:, ['method', 'rank', 'num_phes', 'phes']]

    # print('stats_rank: {}'.format(stats_rank.head()))
    stats_rank = stats_rank.sort_values(['method', 'rank'])
    return stats_rank


def cvs_to_median(stats: pd.DataFrame, stat_col='stat', cols_fix=['method', 'phe']):
    # cols_index = ['method', 'phe']
    cols_index = cols_fix
    # stats_median = stats.groupby(cols_index)['stat'].agg(np.median)
    stats_median = stats.groupby(cols_index)[stat_col].agg(np.median)
    stats_median = stats_median.reset_index(drop=False)
    return stats_median


def run_acc(droot, pformats, dresult, methodms, method_modeld, model_integrate, methodd_type, phes, cvn, fymld):
    datasets = ['ts']

    for max_on in ['nagelkerke']:
        #for plot_on in ['nagelkerke', 'auc', 'auprc']:
        #for plot_on in ['nagelkerke', 'auc']:
        #for plot_on in ['auprc']:
        for plot_on in ['nagelkerke', 'auc', 'auprc']:
            for (withcov, withcov_at_maxon) in get_withcovs(max_on, plot_on):
                logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

                stats = load_best_stats_acc(dresult, methodms, method_modeld, None, None, model_integrate, phes, cvn,
                                            fymld, max_on, plot_on, withcov_at_maxon, withcov, datasets, None)
                logger.debug('stats_acc\n{}'.format(stats.head()))

                # duplicate for plot_on
                dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)

                stats.to_csv(dval / (plot_on + '.tsv'), na_rep='NaN', sep='\t', index=False)

                stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','cvi','stat_on','stat'])
                stats_paper_public.to_csv(dval / (plot_on + '_public.tsv'), na_rep='NaN', sep='\t', index=False)


                # plot
                paper_validate.validate_acc(dval, pformats, stats, plot_on, methodms, methodd_type, phes)

                # stats
                # relative acc tsv
                # stats_covonly = load_best_stats_acc(dresult, ['covonly'], method_modeld, None, None, model_integrate, phes, cvn,
                #                                    fymld, max_on, plot_on, withcov_at_maxon, withcov, datasets)
                # logger.debug('stats_covonly: {}'.format(stats_covonly.head()))

                # TMP: precisely, if boosting methods are 3
                if (plot_on == 'nagelkerke') and len(methodms) == 10:
                    # exclude boosting_nonadd, boosting_add
                    stats_integ = stats.loc[(stats['method'] != 'boosting_add') & (stats['method'] != 'boosting_nonadd'), :]

                    stats_median = cvs_to_median(stats_integ)

                    diff_stats_median = convert_median_relative_stat(stats_median, kind='diff')
                    diff_stats_median = diff_stats_median.reset_index(drop=False)
                    diff_stats_median = diff_stats_median.sort_values(['phe', 'stat'], ascending=[True, False])
                    diff_stats_median.to_csv(dval / (plot_on + '.median_diff_sort.tsv'), na_rep='NaN', sep='\t', index=False)

                    # relave improvement
                    relative_stats_median = convert_median_relative_stat(stats_median)
                    # diff_stats = convert_diff_stat(stats, stats_covonly)
                    relative_stats_median = relative_stats_median.reset_index(drop=False)
                    logger.debug('relative_stats_median: {}'.format(relative_stats_median.head()))
                    relative_stats_median.to_csv(dval / (plot_on + '.median_relative.tsv'), na_rep='NaN', sep='\t', index=False)

                    relative_stats_median = relative_stats_median.sort_values(['phe', 'stat'], ascending=[True, False])
                    relative_stats_median.to_csv(dval / (plot_on + '.median_relative_sort.tsv'), na_rep='NaN', sep='\t', index=False)

                    # stats_top = calculate_top_number(relative_stats_median)
                    # print('stats_top', stats_top)
                    # stats_top.to_csv(dval / (plot_on + '.count_top.tsv'), na_rep='NaN', sep='\t', index=False)

                    stats_rank = calculate_rank(relative_stats_median)
                    # print('stats_rank', stats_rank)
                    stats_rank.to_csv(dval / (plot_on + '.rank.tsv'), na_rep='NaN', sep='\t', index=False)

                    # relative compared to add, or nonadd
                    stats_median = cvs_to_median(stats)

                    for baseline in ['boosting_add', 'boosting_nonadd']:
                        relative_stats_median = convert_median_relative_stat(stats_median, baseline_method=baseline)
                        relative_stats_median = relative_stats_median.reset_index(drop=False)
                        relative_stats_median = relative_stats_median.sort_values(['phe', 'stat'], ascending=[True, False])
                        relative_stats_median.to_csv(dval / (plot_on + '.relative-to-' + baseline +
                                                     '.median_relative_sort.tsv'), na_rep='NaN', sep='\t', index=False)


# TODO: integrate to run_acc
def run_acc_noneur(droot, pformats, dresult, methodms, method_modeld, methodd_mode_ts, model_integrate, methodd_type, phes, cvn, fymld):
    # load allstat from mode=None
    methodd_mode_va = defaultdict(lambda: None)

    datasets = ['ts']

    for max_on in ['nagelkerke']:
        for plot_on in ['nagelkerke', 'auc']:
            for (withcov, withcov_at_maxon) in get_withcovs(max_on, plot_on):
                logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

                stats = load_best_stats_acc(dresult, methodms, method_modeld, methodd_mode_va, methodd_mode_ts, model_integrate, phes, cvn,
                                            fymld, max_on, plot_on, withcov_at_maxon, withcov, datasets, None)
                logger.debug('stats_auc\n{}'.format(stats.head()))

                # duplicate for plot_on
                dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)

                stats.to_csv(dval / (plot_on + '.tsv'), na_rep='NaN', sep='\t', index=False)

                stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','cvi','stat_on','stat'])
                stats_paper_public.to_csv(dval / (plot_on + '_public.tsv'), na_rep='NaN', sep='\t', index=False)


                # plot
                paper_validate.validate_acc(dval, pformats, stats, plot_on, methodms, methodd_type, phes)
                # paper_validate.validate_acc(dval, pformats, stats, plot_on, methodms, phes)
                # validate_acc(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, plot_on, withcov_at_maxon, withcov)


# TODO: integrate to run_acc
def run_acc_noneur_covva(droot, pformats, dresult, methodms, method_modeld, methodd_mode_ts, model_integrate, methodd_type, phes, cvn, fymld):
    # load allstat from mode=None
    methodd_mode_va = defaultdict(lambda: None)

    # TMP
    # select same parameter as eur
    method_modeld_va = method_modeld.copy()
    for method in method_modeld_va:
        if method in 'lassosum_noneur,ldpred_noneur,prscs_noneur,sbayesr_noneur,clump_noneur'.split(','):
            method_cor = method.split('_')[0]
            method_modeld_va[method] = method_modeld[method_cor]

    datasets = ['ts']

    for max_on in ['nagelkerke']:
        for plot_on in ['nagelkerke', 'auc']:
            for (withcov, withcov_at_maxon) in get_withcovs(max_on, plot_on):
                logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

                stats = load_best_stats_acc(dresult, methodms, method_modeld, methodd_mode_va, methodd_mode_ts, model_integrate, phes, cvn,
                                            fymld, max_on, plot_on, withcov_at_maxon, withcov, datasets, method_modeld_va)
                logger.debug('stats_auc\n{}'.format(stats.head()))

                # duplicate for plot_on
                dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)

                stats.to_csv(dval / (plot_on + '.tsv'), na_rep='NaN', sep='\t', index=False)

                stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','cvi','stat_on','stat'])
                stats_paper_public.to_csv(dval / (plot_on + '_public.tsv'), na_rep='NaN', sep='\t', index=False)


                # plot
                paper_validate.validate_acc(dval, pformats, stats, plot_on, methodms, methodd_type, phes)


# TODO: integrate to run_acc
def run_acc_noneur_covva_boost(droot, pformats, dresult, methodms, method_modeld, methodd_mode_ts, model_integrate, methodd_type, phes, cvn, fymld):
    # load allstat from mode=None
    methodd_mode_va = defaultdict(lambda: None)

    # TMP
    # make all
    # select same parameter as eur
    method_modeld_va = method_modeld.copy()
    # for method in method_modeld_va:
    #    if method in 'lassosum_noneur,ldpred_noneur,prscs_noneur,sbayesr_noneur,clump_noneur'.split(','):
    #        method_cor = method.split('_')[0]
    #        method_modeld_va[method] = method_modeld[method_cor]

    datasets = ['ts']

    for max_on in ['nagelkerke']:
        for plot_on in ['nagelkerke', 'auc']:
            for (withcov, withcov_at_maxon) in get_withcovs(max_on, plot_on):
                logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

                stats = load_best_stats_acc(dresult, methodms, method_modeld, methodd_mode_va, methodd_mode_ts, model_integrate, phes, cvn,
                                            fymld, max_on, plot_on, withcov_at_maxon, withcov, datasets, method_modeld_va)
                logger.debug('stats_auc\n{}'.format(stats.head()))

                # duplicate for plot_on
                dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)

                stats.to_csv(dval / (plot_on + '.tsv'), na_rep='NaN', sep='\t', index=False)

                stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','cvi','stat_on','stat'])
                stats_paper_public.to_csv(dval / (plot_on + '_public.tsv'), na_rep='NaN', sep='\t', index=False)


                # plot
                paper_validate.validate_acc(dval, pformats, stats, plot_on, methodms, methodd_type, phes)

                methods_noneur = [x for x in methodms if 'noneur' in x]
                stats = stats.loc[stats['method'].isin(methods_noneur), :]

                if plot_on == 'nagelkerke':
                    # def cvs_to_median(stats: pd.DataFrame):
                    #    cols_index = ['method', 'phe']
                    #    stats_median = stats.groupby(cols_index)['stat'].agg(np.median)
                    #    stats_median = stats_median.reset_index(drop=False)
                    #    return stats_median
                    # exclude boosting_nonadd, boosting_add
                    stats_integ = stats.loc[(stats['method'] != 'boosting_add-noneur') & (stats['method'] != 'boosting_nonadd-noneur'), :]

                    stats_median = cvs_to_median(stats_integ)
                    logger.debug('stats_median: {}'.format(stats_median))
                    relative_stats_median = convert_median_relative_stat(stats_median, baseline_method='boosting_integrate-noneur')
                    # relative_stats = convert_relative_stat(stats, stats_covonly)
                    relative_stats_median = relative_stats_median.reset_index(drop=False)
                    logger.debug('relative_stats_median: {}'.format(relative_stats_median.head()))

                    relative_stats_median.to_csv(dval / (plot_on + '.median_relative.tsv'), na_rep='NaN', sep='\t', index=False)

                    relative_stats_median = relative_stats_median.sort_values(['phe', 'stat'], ascending=[True, False])
                    relative_stats_median.to_csv(dval / (plot_on + '.median_relative_sort.tsv'), na_rep='NaN', sep='\t', index=False)

                    # stats_top = calculate_top_number(relative_stats_median)
                    # print('stats_top', stats_top)
                    # stats_top.to_csv(dval / (plot_on + '.count_top.tsv'), na_rep='NaN', sep='\t', index=False)

                    stats_rank = calculate_rank(relative_stats_median)
                    # print('stats_rank', stats_rank)
                    stats_rank.to_csv(dval / (plot_on + '.rank.tsv'), na_rep='NaN', sep='\t', index=False)

                    # raise NotImplementedError


def run_acc_simu(droot, pformats, ddata, dresult, methodms, method_modeld, model_integrate, simu_types,
                 cvn, fymld):
    datasets = ['ts']

    for max_on in ['nagelkerke']:
        for plot_on in ['nagelkerke', 'auc', 'auprc']:
            # for plot_on in ['nagelkerke', 'auc']:
            # exception
            withcovs = [(False, False)]
            for (withcov, withcov_at_maxon) in withcovs:
                # for (withcov, withcov_at_maxon) in get_withcovs(max_on, plot_on):
                logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

                stats = load_best_stats_acc_simu(ddata, dresult, methodms, method_modeld, model_integrate, simu_types,
                                                 cvn,
                                                 fymld, max_on, plot_on, withcov_at_maxon, withcov, datasets)
                logger.debug('stats_acc\n{}'.format(stats.head()))

                # duplicate for plot_on
                dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)
                # TODO: write down summary stats
                # save to dval / 'summary' / .tsv

                stats.to_csv(dval / (plot_on + '.tsv'), na_rep='NaN', sep='\t', index=False)

                stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','h2add','h2dom','ncausal','stat_on','stat'])
                stats_paper_public.to_csv(dval / (plot_on + '_public.tsv'), na_rep='NaN', sep='\t', index=False)


                # plot
                paper_validate.validate_acc_simu(dval, pformats, stats, plot_on, methodms)
                # validate_acc(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, plot_on, withcov_at_maxon, withcov)


def run_odds_ratio(droot, pformats, dresult, methodms, method_modeld, model_integrate, fpheno, phes, phes_female, fcv, cvn, fymld):
    or_props = [1.0, 3.0, 5.0, 10.0]
    # or_props = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0]
    datasets = ['ts']

    for max_on in ['nagelkerke']:
        # withcovs_and_maxon = get_withcovs(max_on)
        for (withcov, withcov_at_maxon) in get_withcovs(max_on):
            logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))
            stats = load_best_or(dresult, methodms, method_modeld, model_integrate,
                                 fpheno, phes, phes_female, fcv, cvn,
                                 fymld, max_on, withcov_at_maxon, withcov,
                                 datasets, or_props)
            logger.debug('stats\n{}'.format(stats.head()))
            logger.debug('stats\n{}'.format(stats.columns))

            # plot
            dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)

            stats.to_csv(dval / ( 'stat.tsv'), na_rep='NaN', sep='\t', index=False)

            stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','cvi','stat_on','stat'])
            stats_paper_public.to_csv(dval / ( 'stat_public.tsv'), na_rep='NaN', sep='\t', index=False)
        

            paper_validate.validate_odds_ratio(dval, pformats, stats, methodms, phes)

            stats_median = cvs_to_median(stats, cols_fix=['method', 'phe', 'stat_on'])
            stats_median = stats_median.sort_values(['phe', 'stat_on', 'stat'], ascending=[True, True, False])
            stats_median.to_csv(dval / ('or.median.tsv'), na_rep='NaN', sep='\t', index=False)


def get_phes_bothsex(phes, phes_female):
    return [x for x in phes if x not in phes_female]


def get_phes_female(phes, phes_female):
    return [x for x in phes if x in phes_female]


def load_phes_correlation(
        fpheno, phes, phes_female,
        sex, fcv, cvi,
        dataset, stat_name, stat_fn):

    cv = sampleio.load_cv(fcv)
    pheno = sampleio.load_pheno(fpheno)

    if sex == 'both':
        phes_sex = get_phes_bothsex(phes, phes_female)
    elif sex == 'female':
        # use bothsex+female
        phes_sex = phes
        # phes_sex = get_phes_female(phes, phes_female)
    else:
        raise NotImplementedError

    pheno = sampleop.extract_dataset(pheno, dataset=dataset, cvi=cvi, cvs=cv, sex=sex)
    # pheno = sampleop.extract_dataset(pheno, dataset='ts', cvi=cvi, cvs=cv, sex=sex)

    pheno_use = pheno.loc[:, phes_sex].copy()
    # pheno_use = pheno.loc[samples, phes_sex].copy()

    stats = pd.DataFrame([], index=phes_sex, columns=phes_sex, dtype=float)
    for phe1 in phes_sex:
        for phe2 in phes_sex:
            stat = stat_fn(pheno_use[phe1], pheno_use[phe2])
            # stat = jaccard_similarity(pheno_use[phe1], pheno_use[phe2])
            stats.loc[phe1, phe2] = stat

    return stats


def jaccard_similarity(phes1, phes2):
    stat = metrics.jaccard_score(phes1, phes2)
    return stat


def chisq_test(phes1, phes2):
    res = scistats.contingency.crosstab(phes1, phes2)
    count = pd.DataFrame(res.count, index=res.elements[0], columns=res.elements[1])
    # format order
    label_order = [0, 1]
    count = count.loc[label_order, label_order]
    # logger.debug('count: {}'.format(count))

    stat = scistats.chi2_contingency(count)
    stat = stat.pvalue
    # logger.debug('stat: {:e}'.format(stat))

    # do this in plot
    # if stat==0.0:
    # cannot use np.float128, so should use smallest float
    # TMP
    # stat= np.nextafter(0,1)

    return stat


def run_phes_correlation(droot, pformats, dresult, methodms, method_modeld, model_integrate, fpheno, phes, phes_female, fcv, cvis, fymld):
    dataset = 'tr'
    # dataset='ts'

    stat_fns = [('jaccard_similarity', jaccard_similarity),
                ('chisq', chisq_test)]

    for max_on in ['nagelkerke']:
        for (withcov, withcov_at_maxon) in get_withcovs(max_on, nocov_only=True):
            logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

            for cvi in cvis:
                for sex in ['both', 'female']:
                    for (stat_name, stat_fn) in stat_fns:

                        stat = load_phes_correlation(fpheno, phes, phes_female, sex, fcv, cvi, dataset, stat_name, stat_fn)

                        # plot
                        # dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)
                        dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon, cvi=cvi)
                        paper_validate.validate_phes_correlation(dval, pformats, stat, sex, stat_name)
                        # validate_paper.auc.validate_auc(droot, stats)


def load_h2(dresult, phes, phes_female, cvi=None, kind='default'):

    # stats = pd.DataFrame([], index=phes, columns=phes_sex, dtype=float)
    stats = []
    for phe in phes:
        if kind == 'default':
            fh2 = dresult / 'd-ldsc' / kind / phe / ('cv' + str(cvi)) / 'd-ldsc.h2'
        elif kind == 'all-sample' or kind == 'neale-lab' or kind == 'neale-lab_all-snvs' or kind=='neale-lab_max-chisq':
            fh2 = dresult / 'd-ldsc' / kind / phe / 'd-ldsc.h2'
        else:
            raise RuntimeError('Unknown kind: ',kind)

        if not os.path.exists(fh2):
            logger.debug('not exists: {}'.format(fh2))
            continue
        stat = pd.read_csv(fh2, sep='\t', header=0)
        logger.debug('stat\n{}'.format(stat))
        stat.loc[:, 'phe'] = phe
        stats += [stat]
    stats = pd.concat(stats, axis=0, ignore_index=True)

    logger.debug('stats\n{}'.format(stats.all()))
    logger.debug('stats\n{}'.format(stats.columns))
    assert len(stats.columns) == 29

    stats.loc[:, 'prop_h2'] = stats['h2_liability_d'] / stats['h2_liability']

    # h2_liability
    # logger.debug('stats\n{}'.format(stats.loc[:, ['phe', 'h2_observed', 'h2_liability', 'h2_observed_d', 'h2_liability_d']]))
    # logger.debug('stats\n{}'.format(stats.loc[:, ['phe', 'h2_observed', 'h2_observed_se',
    #             'h2_observed_d', 'h2_observed_se_d', 'h2_p', 'h2_p_d']].to_string()))
    logger.debug('stats\n{}'.format(stats.loc[:, ['phe', 'prop_h2', 'h2_liability', 'h2_liability_se',
                 'h2_liability_d', 'h2_liability_se_d', 'h2_p', 'h2_p_d']].to_string()))

    # for phe1 in phes_sex:
    #    for phe2 in phes_sex:
    #        stat = stat_fn(pheno_use[phe1], pheno_use[phe2])
    #        # stat = jaccard_similarity(pheno_use[phe1], pheno_use[phe2])
    #        stats.loc[phe1, phe2] = stat

    # raise NotImplementedError

    return stats


# cv
def run_h2(droot, pformats, dresult, fpheno, phes, phes_female, fcv, cvis, fymld):
    dataset = 'tr'
    # dataset='ts'

    for max_on in ['nagelkerke']:
        for (withcov, withcov_at_maxon) in get_withcovs(max_on, nocov_only=True):
            logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

            for cvi in cvis:
                # for sex in ['both', 'female']:
                # for (stat_name, stat_fn) in stat_fns:

                stat = load_h2(dresult, phes, phes_female, cvi)
                # stat = load_h2(dresult, fpheno, phes, phes_female, fcv, cvi, dataset)

                # plot
                # dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)
                dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon, cvi=cvi)
                paper_validate.validate_h2(dval, pformats, stat, sex, stat_name)
                # validate_paper.auc.validate_auc(droot, stats)


# without cv
def run_h2_all_sample(droot, pformats, dresult, fpheno, phes, phes_female, fcv, cvis, fymld, kind='all-sample'):
    dataset = 'tr'
    # dataset='ts'

    for max_on in ['nagelkerke']:
        for (withcov, withcov_at_maxon) in get_withcovs(max_on, nocov_only=True):
            logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))

            # for cvi in cvis:
            # for sex in ['both', 'female']:
            # for (stat_name, stat_fn) in stat_fns:

            stat = load_h2(dresult, phes, phes_female, kind=kind)
            # stat = load_h2(dresult, phes, phes_female, kind='all-sample')

            # plot
            dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)
            # dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon, cvi=cvi)
            paper_validate.validate_h2(dval, pformats, stat)
            # validate_paper.auc.validate_auc(droot, stats)


def run_sample_score(droot, pformats, dresult, methodms, method_modeld, model_integrate, fpheno, phes, phes_female, fcv, cvis, fymld):
    datasets = ['ts']

    for max_on in ['nagelkerke']:
        # withcovs_and_maxon = get_withcovs(max_on)
        # for (withcov, withcov_at_maxon) in get_withcovs(max_on):
        for (withcov, withcov_at_maxon) in get_withcovs(max_on, nocov_only=True):
            # TMP: do in get_withcovs
            # if withcov:
            # continue
            logger.debug('withcov, withcov_at_maxon: {}, {}'.format(withcov, withcov_at_maxon))
            for phe in phes:
                sex = sex_of_phe(phe, phes_female)
                for cvi in cvis:
                    scored = load_sample_scored(dresult, methodms, method_modeld, model_integrate,
                                                fpheno, phe, sex, fcv, cvi,
                                                fymld, max_on, withcov_at_maxon, withcov, datasets)

                    # plot
                    # dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon)
                    dval = dname_maxon_create(droot, max_on, withcov, withcov_at_maxon, phe, cvi)
                    paper_validate.validate_sample_score(dval, pformats, scored, methodms, phe)
                    # validate_paper.auc.validate_auc(droot, stats)


def run_nsnv(droot, pformats, dresult, methodms, method_modeld, model_integrate, phes, cvn, fymld):

    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))
            stats = load_best_nsnv(dresult, methodms, method_modeld, None, model_integrate, phes, cvn, fymld, max_on, withcov_at_maxon)
            logger.debug('stats_nsnv\n{}'.format(stats.head()))

            # plot
            dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon)

            stats.to_csv(dval / ( 'stat.tsv'), na_rep='NaN', sep='\t', index=False)

            stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['method','phe','cvi','nsnv_use'])
            stats_paper_public.to_csv(dval / ( 'stat_public.tsv'), na_rep='NaN', sep='\t', index=False)
            

            paper_validate.validate_nsnv(dval, pformats, stats, methodms, phes)

            stats_median = cvs_to_median(stats, stat_col='nsnv_use')
            stats_median = stats_median.sort_values(['phe', 'nsnv_use'], ascending=[True, True])
            stats_median.to_csv(dval / ('nsnv.median.tsv'), na_rep='NaN', sep='\t', index=False)


def run_genetic_model(droot, pformats, dresult, method_modeld, phes, cvis, fymld):

    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))
            for cvi in cvis:
                stats = load_genetic_model(dresult, method_modeld, phes, cvi, fymld, max_on, withcov_at_maxon)
                logger.debug('stats\n{}'.format(stats.head()))

                



                # plot
                dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon, cvi=cvi)

                stats.to_csv(dval / ('stat.tsv'), na_rep='NaN', sep='\t', index=False)
                stats_paper_public=paper_public_stat.stat_to_paper_public(stats,['phe','cvi','nsnv_use','overdom','dom','add','rec','overrec'])
                stats_paper_public.to_csv(dval / ( 'stat_public.tsv'), na_rep='NaN', sep='\t', index=False)

                paper_validate.validate_genetic_model(dval, pformats, stats, phes)
                # validate_paper.auc.validate_auc(droot, stats)


                # validate_nsnv(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, withcov_at_maxon)


def run_genetic_model_gwas(droot, pformats, dresult, method_modeld, phes, cvis, fymld):

    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))
            for phe in phes:
                for cvi in cvis:
                    stat = load_genetic_model_gwas(dresult, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon)
                    # logger.debug('stats\n{}'.format(stats.head()))

                    # plot
                    dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon, phe=phe, cvi=cvi)

                    # in validate.py
                    #stat.to_csv(dval / ('stat.tsv'), na_rep='NaN', sep='\t', index=False)
                    #stat_paper_public=paper_public_stat.stat_to_paper_public(stat,['phe','cvi','nsnv_use','overdom','dom','add','rec','overrec'])
                    #stat_paper_public.to_csv(dval / ( 'stat_public.tsv'), na_rep='NaN', sep='\t', index=False)
                    

                    paper_validate.validate_genetic_model_gwas(dval, pformats, stat, phe)
                    # validate_paper.auc.validate_auc(droot, stats)

                    # validate_nsnv(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, withcov_at_maxon)


def run_genetic_model_vs_nsnv(droot, pformats, dresult, method_modeld, phes, cvis, fymld):

    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))
            for cvi in cvis:
                # genetic_model already has 'nsnv_use' in stats
                stats = load_genetic_model(dresult, method_modeld, phes, cvi, fymld, max_on, withcov_at_maxon)
                logger.debug('stats\n{}'.format(stats.head()))

                # plot
                dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon, cvi=cvi)
                paper_validate.validate_genetic_model_vs_nsnv(dval, pformats, stats, phes)
                # validate_paper.auc.validate_auc(droot, stats)

                # validate_nsnv(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, withcov_at_maxon)


def load_nonaddlist(dresult, methodm, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon, gmodels):
    '''
    Load ss containing loss, gwas of selected snvs in wgt
    '''
    regond = load_regon(fymld)
    # stats = []

    # ssd = load_ssd_raw(dresult, phe, cvi, ['add'])
    # ssd = ssop.load_ssd_raw(dresult, phe, cvi, gmodels)

    # methodm = 'boosting_nonadd'
    reg_on, va_on = regon_vaon(methodm, regond, withcov_at_maxon)
    paras_best_cvs = load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                        va_on, fymld['cand'], method_modeld, None, cvis=[cvi])

    wgt = load_best_wgt(dresult, methodm, phe, cvi, paras_best_cvs, method_modeld)
    ssd = update_ssd_wgt(ssd, wgt, methodm, use_sid_tmp=('dataset.12' in str(dresult)))

    # for loss
    # methodm = 'boosting_loss'
    # wgt = load_best_wgt(dresult, methodm, phe, cvi, None, method_modeld)
    # ssd = update_ssd_wgt(ssd, wgt, methodm, use_sid_tmp=('dataset.12' in str(dresult)))

    # stats.extend(stats_best_cvs)

    # remove None
    # stats = [x for x in stats if x is not None]
    # stats = pd.DataFrame.from_records(stats)
    return ssd


def run_manhattan(droot, pformats, dresult, methodms, method_modeld, phes, cvis, fymld):
    # logger.debug('cvn {}'.format(cvn))

    gmodels = ['add']

    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))
            for phe in phes:
                # for cvi in range(cvn):
                for cvi in cvis:
                    ssd, methodms_use = load_manhattan_ssd(dresult, methodms, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon, gmodels)
                    # logger.debug('stats_auc\n{}'.format(ssd.head()))
                    logger.debug('methodms_use: {}'.format(methodms_use))

                    # plot
                    dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon, phe=phe, cvi=cvi)
                    paper_validate.validate_manhattan(dval, pformats, ssd, methodms_use, phe, gmodels=['add'])
                    # paper_validate.validate_manhattan(dval, pformats, ssd, methodms, phe, gmodels=['add'])
                    # validate_paper.auc.validate_auc(droot, stats)

            # validate_nsnv(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, withcov_at_maxon)


def locus_zoom_ranges():
    local_rangesd = {
        'psr': [
            #[
            #    ((22, 20100000, 21500000), ('dom', 'add'), ['22:20876509:G:A']),
            #    ((2, 110800000, 112200000), ('rec', 'add'), ['2:111600519:T:C']),
            #    ((6, 167700000, 169100000), ('hetonly', 'add'), ['6:168413393:A:G']),
            #    ((2, 228200000, 229600000), ('rec', 'add'), ['2:228968345:G:A']),
            #    ((2, 64800000, 66200000), ('hetonly', 'add'), ['2:65586769:C:A']),
            #    ], 
                [
                    ((11, 95800000, 99200000), ('add', 'rec', 'hetonly'), [
                     '11:96523529:T:C', '11:97281245:A:G', '11:97734972:C:T', '11:98176280:C:T', ], False),  # for loss; 11:96523529
        ]],  # for loss; 11:96523529
        'ra': [[
            ((18, 59300000, 60700000), ('dom', 'add'), ['18:60004253:G:A']),
            ((2, 143600000, 145000000), ('dom', 'add'), ['2:144385104:C:T']),
            ((3, 118400000, 119900000), ('hetonly', 'add'), ['3:119123814:C:T']),
        ]],
    }

    return local_rangesd


def run_locus_zoom(droot, pformats, dresult, methodms, method_modeld, phes, cvis, fymld, fgene_exon):
    gmodels = ['add', 'dom', 'rec', 'hetonly']
    gene_exon = pd.read_csv(fgene_exon, sep='\t')
    logger.debug('gene_exon: {}'.format(gene_exon.head()))
    # raise NotImplementedError

    local_rangesd = locus_zoom_ranges()
    for max_on in ['nagelkerke']:
        for withcov_at_maxon in get_withcov_at_maxon(max_on):
            logger.debug('withcov_at_maxon: {}'.format(withcov_at_maxon))
            for phe in phes:
                if phe not in local_rangesd:
                    continue
                for cvi in cvis:
                    ssd, _ = load_manhattan_ssd(dresult, methodms, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon, gmodels)

                    # plot
                    dval = dname_maxon_create(droot, max_on, withcov_at_maxon=withcov_at_maxon, phe=phe, cvi=cvi)
                    paper_validate.validate_locus_zoom(dval, pformats, ssd, methodms, phe, local_rangesd[phe], gene_exon)
                    # validate_paper.auc.validate_auc(droot, stats)

            # validate_nsnv(dval, pformats, dresult, methodms, method_modeld, phes, cvn, fymld, max_on, withcov_at_maxon)


def extract_boosing_methods(methods, ext):
    methods_ext = []
    for method in methods:
        if method.startswith('boosting_') and method not in ext:
            continue
        else:
            methods_ext.append(method)
    return methods_ext


def paper(dout, pformats, runs, dresult, ddata,
          methodms, methodd_model, methodd_mode, methodd_integrate, methodd_type,
          fpheno, phes, phes_female, fcv, cvn,
          para_cand, para_regon, para_com,
          fgene_exon):

    logger.debug('runs: {}'.format(runs))

    fymld = get_fymld(ddata, para_cand, para_regon, para_com)

    # TODO: runs = ['best-para','nonadd-list',
    # 'compare-gwas','genetic-model-aic',
    # 'auc-split']

    # TODO: make runs in bash
    # runs = ['genetic-model-gwas', 'nsnv']

    # DONE
    # runs = ['sample-score']
    # runs = ['manhattan']
    # runs = ['locus-zoom']
    # runs = ['genetic-model']
    # runs = ['nsnv']
    # runs = ['acc']
    # runs = ['odds-ratio']

    # runs = ['sample-score', 'manhattan', 'locus-zoom', 'genetic-model', 'nsnv', 'acc', 'odds-ratio']
    # runs = ['genetic-model', 'nsnv', 'acc', 'odds-ratio']
    # runs = ['acc', 'odds-ratio', 'nsnv']
    # runs = ['acc']
    # runs = ['acc-simu']
    # runs = ['acc-noneur']
    # runs = ['para-grid']
    # runs = ['genetic-model-vs-nsnv']
    # runs = ['phes-correlation']

    if 'best-para' in runs:
        logger.debug('run best-para')
        droot = dout / 'best-para'
        run_best_para_table(droot, pformats, dresult, methodms, methodd_model, methodd_integrate, phes, cvn, fymld)

    if 'para-grid' in runs:
        logger.debug('run para-grid')
        droot = dout / 'para-grid'
        methodms_ext = ['boosting_nonadd', 'boosting_add']
        run_para_grid(droot, pformats, dresult, methodms_ext, methodd_model, methodd_integrate, phes, [0], fymld)

    if 'acc' in runs:
        logger.debug('run acc')
        droot = dout / 'acc'
        run_acc(droot, pformats, dresult, methodms, methodd_model, methodd_integrate, methodd_type, phes, cvn, fymld)

    if 'acc-simu' in runs:
        logger.debug('run acc-simu')

        # use .txt?
        simu_types = ["h2add-0.2_h2dom-0.05_ncausal-100",
                      "h2add-0.2_h2dom-0.05_ncausal-1000",
                      "h2add-0.25_h2dom-0.0_ncausal-100",
                      "h2add-0.25_h2dom-0.0_ncausal-1000",
                      "h2add-0.08_h2dom-0.02_ncausal-100",
                      "h2add-0.08_h2dom-0.02_ncausal-1000",
                      "h2add-0.1_h2dom-0.0_ncausal-100",
                      "h2add-0.1_h2dom-0.0_ncausal-1000"]

        droot = dout / 'acc-simu'
        run_acc_simu(droot, pformats, ddata, dresult, methodms, methodd_model, methodd_integrate, simu_types, 1, fymld)
        # run_acc_simu(droot, pformats, ddata, dresult, methodms, methodd_model, methodd_integrate, simu_types, phes, 1, fymld)

    if 'acc-split-dom' in runs:
        logger.debug('run acc-split-dom')
        droot = dout / 'acc-split-dom'
        # run_acc_split_dom(droot, pformats, dresult, methodms, methodd_model, phes, cvn, fymld)

    # old
    #if 'acc-noneur' in runs:
    #    logger.debug('run acc-noneur')
    #    droot = dout / 'acc-noneur'
    #    run_acc_noneur(droot, pformats, dresult, methodms, methodd_model, methodd_mode, methodd_integrate, methodd_type, phes, cvn, fymld)

    # old
    #if 'acc-noneur-covva' in runs:
    #    # TMP
    #    logger.debug('run acc-noneur_covva')
    #    droot = dout / 'acc-noneu-covvar'
    #    run_acc_noneur_covva(droot, pformats, dresult, methodms, methodd_model, methodd_mode, methodd_integrate, methodd_type, phes, cvn, fymld)

    if 'acc-noneur-covva-boost' in runs:
        # TMP
        logger.debug('run acc-noneur_covva-boost')
        droot = dout / 'acc-noneu-covvar-boost'
        run_acc_noneur_covva_boost(droot, pformats, dresult, methodms, methodd_model, methodd_mode, methodd_integrate, methodd_type, phes, cvn, fymld)

    if 'odds-ratio' in runs:
        logger.debug('run odds ratio')
        droot = dout / 'odds-ratio'
        run_odds_ratio(droot, pformats, dresult, methodms, methodd_model, methodd_integrate, fpheno, phes, phes_female, fcv, cvn, fymld)

    if 'nsnv' in runs:
        # double count variants if variant names are different (even when chrompos is the same)
        logger.debug('run nsnv')
        droot = dout / 'nsnv'
        run_nsnv(droot, pformats, dresult, methodms, methodd_model, methodd_integrate, phes, cvn, fymld)

    if 'genetic-model' in runs:
        logger.debug('run genetic model')
        droot = dout / 'genetic-model'
        if 'boosting_nonadd' in methodms:
            run_genetic_model(droot, pformats, dresult, methodd_model, phes, [0], fymld)

    if 'genetic-model-gwas' in runs:
        logger.debug('run genetic model gwas')
        droot = dout / 'genetic-model-gwas'
        # TMP
        #phes_ext = ['atm', 'ra', 'psr']
        phes_ext = [ 'psr']
        if 'boosting_nonadd' in methodms:
            run_genetic_model_gwas(droot, pformats, dresult, methodd_model, phes_ext, [0], fymld)

    if 'genetic-model-vs-nsnv' in runs:
        logger.debug('run genetic model vs nsnv')
        droot = dout / 'genetic-model-vs-nsnv'
        if 'boosting_nonadd' in methodms:
            run_genetic_model_vs_nsnv(droot, pformats, dresult, methodd_model, phes, [0], fymld)

    # in analysis.py
    # if 'nonaddlist' in runs:
    #    logger.debug('run nonaddlist')
    #    phes_ext = [ 'ra', 'psr']
    #    droot = dout / 'nonaddlist'
    #    if 'boosting_nonadd' in methodms:
    #        methodm_ext ='boosting_nonadd'
    #        run_nonaddlist(droot, pformats, dresult, methodm, methodd_model, phes_ext, [0], fymld)

    if 'manhattan' in runs:
        logger.debug('run manhattan')
        phes_ext = ['atm', 'ra', 'psr']
        # phes_ext = ['atm']
        droot = dout / 'manhattan'
        if 'boosting_nonadd' in methodms:
            # 'boosting_nonadd' only here
            methodms_ext = extract_boosing_methods(methodms, ['boosting_nonadd', 'boosting_add'])
            # TMP
            # methodms_ext = 'boosting_nonadd,ldpred'.split(',')
            run_manhattan(droot, pformats, dresult, methodms_ext, methodd_model, phes_ext, [0], fymld)

    if 'locus-zoom' in runs:
        logger.debug('run locus zoom')
        phes_ext = ['psr']
        # phes_ext = ['ra', 'psr']
        droot = dout / 'locus-zoom'
        if 'boosting_nonadd' in methodms:
            # 'boosting_nonadd' only here
            methodms_ext = 'boosting_nonadd'.split(',')
            run_locus_zoom(droot, pformats, dresult, methodms_ext, methodd_model, phes_ext, [0], fymld, fgene_exon)
            # run_genetic_model(droot, pformats, dresult, method_modeld, phes, cvn, fymld)

    if 'phes-correlation' in runs:
        logger.debug('run phes correlation')
        droot = dout / 'phes-correlation'
        run_phes_correlation(droot, pformats, dresult, methodms, methodd_model, methodd_integrate, fpheno, phes, phes_female, fcv, [0], fymld)

    if 'h2' in runs:
        logger.debug('run h2')
        droot = dout / 'h2'
        run_h2(droot, pformats, dresult, fpheno, phes, phes_female, fcv, [0], fymld)

    if 'h2_all-sample' in runs:
        logger.debug('run h2_all-sample')
        droot = dout / 'h2_all-sample'
        run_h2_all_sample(droot, pformats, dresult, fpheno, phes, phes_female, fcv, [0], fymld)

    if 'h2_neale-lab' in runs:
        logger.debug('run h2_neale-lab')
        droot = dout / 'h2_neale-lab'
        run_h2_all_sample(droot, pformats, dresult, fpheno, phes, phes_female, fcv, [0], fymld, kind='neale-lab')


    if 'h2_neale-lab_max-chisq' in runs:
        logger.debug('run h2_neale-lab_max-chisq')
        droot = dout / 'h2_neale-lab_max-chisq'
        run_h2_all_sample(droot, pformats, dresult, fpheno, phes, phes_female, fcv, [0], fymld, kind='neale-lab_max-chisq')

    if 'h2_neale-lab_all-snvs' in runs:
        logger.debug('run h2_neale-lab')
        droot = dout / 'h2_neale-lab_all-snvs'
        run_h2_all_sample(droot, pformats, dresult, fpheno, phes, phes_female, fcv, [0], fymld, kind='neale-lab_all-snvs')

    if 'sample-score' in runs:
        logger.debug('run sample score')
        phes_ext = ['atm']
        # phes_ext = ['atm', 'ra', 'psr']
        droot = dout / 'sample-score'
        run_sample_score(droot, pformats, dresult, methodms, methodd_model, methodd_integrate, fpheno, phes_ext, phes_female, fcv, [0], fymld)
        # run_odds_ratio(droot, pformats, dresult, methodms, method_modeld, fpheno, phes, phes_female, fcv, cvn, fymld)


def main():

    args = argument()

    if args.format is None:
        # pformats=['paper']
        pformats = ['paper', 'slide', 'poster']
    else:
        pformats = args.format

    # 'boosting_integrate' not 'boosting'
    # FIXME: use argparse
    # methodms = args.methods.split(' ')
    # phes = args.phes.split(' ')
    # logger.debug('phes: {}'.format(phes))
    # phes_female = args.phes_female.split(' ')
    # if phes_female is not in phes, then they will not be used
    # if not set(phes_female).issubset(set(phes)):
    #    raise RuntimeError('Wrong phes_female', phes_female)

    # 'boosting_nonadd' and 'boosting_add' are necessary for 'boosting

    # method-> mid_path
    # ex. (kind), (kind, model), (kind, simu_type, model)
    methodd_mid_path = ioarg.arg_method_mid_path(args.method_mid_path)
    logger.debug('method_modeld {}'.format(methodd_mid_path))

    # FIXME: touse
    methodd_type = ioarg.arg_method_type(args.methods_type)
    # logger.debug('method_type {}'.format(methodd_type))

    # method-> (kind, model)
    # method_modeld = arg_method_model_kind(args.method_model)
    # method_modeld = defaultdict(lambda: ('default', None))
    # method_modeld = defaultdict(lambda: None)
    # print('method_model',args.method_model)
    # method_modeld |= dict(args.method_model)
    # method_modeld |= {
    # 'boosting_nonadd': args.model_boost,
    # 'boosting_add': args.model_boost_add,
    # 'boosting_loss': args.model_boost_loss
    # }
    # logger.debug('method_modeld {}'.format(method_modeld))

    # raise NotImplementedError

    methodd_integrate = {
        'boosting_integrate': ['boosting_add', 'boosting_nonadd'],
    }
    if args.model_integrate is not None:
        d = dict(args.model_integrate)
        for k in d:
            d[k] = d[k].split(' ')
        methodd_integrate |= d

    logger.debug('model_integrate {}'.format(methodd_integrate))
    # raise NotImplementedError

    # unnecessary? since we know the noneur is common
    # then, methodm is unknown; or automatically produce? method+'-noneur
    methodd_mode = defaultdict(lambda: None)
    if args.method_mode is not None:
        methodd_mode |= dict(args.method_mode)

    paper(pathlib.Path(args.dout), pformats, args.runs,
          pathlib.Path(args.dresult), pathlib.Path(args.ddata),
          args.methods, methodd_mid_path, methodd_mode,
          methodd_integrate,
          methodd_type,
          pathlib.Path(args.fpheno), args.phes, args.phes_female,
          pathlib.Path(args.fcv), int(args.cvn),
          args.para_cand, args.para_regon, args.para_com,
          args.fgene_exon)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
