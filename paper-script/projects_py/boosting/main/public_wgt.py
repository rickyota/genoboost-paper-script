
import argparse
import os
import pathlib
import yaml
import functools
import pandas as pd
from collections import defaultdict
from logging import getLogger

from ...system.logger import logger_setting
# from ...genetics.sample import sampleio as sampleio
# from ...genetics.score import scoreio as scoreio
from ..io import filename as iof
# from ..score import score as scoreop
from . import paper

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--dout',
                        help='output dir')

    # ./result/nonaddpgs/dataset.12/
    parser.add_argument('--dresult',
                        help='result dir')

    parser.add_argument('--ddata',
                        help='data dir')

    # assume  --method-model $(method) $(model)
    parser.add_argument('--method-model',
                        nargs="+",
                        action='append',
                        help='')

    parser.add_argument('--methods',
                        help='')

    parser.add_argument('--cvn',
                        help='')

    parser.add_argument('--phes',
                        help='')

    parser.add_argument('--para-cand',
                        help='')

    parser.add_argument('--para-com',
                        help='')

    parser.add_argument('--para-regon',
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


# allstat is string
# when withcov=False, regon is not necessary
@functools.cache
def load_allstat(dresult, method, phe, withcov, regon=None, boost_model=None):
    fallstat = iof.file_allstat(dresult, method, phe, withcov, regon, model=boost_model)
    allstat = pd.read_csv(fallstat, sep='\t', index_col=0)
    # para should be string, otherwise, cannot access to filename
    # -> but troublesome for cvi, nsnv etc.
    # allstat = pd.read_csv(fallstat, sep='\t', index_col=0, dtype=str)

    # print('allstat',allstat)
    return allstat


def load_allstatd_integrate(dresult, method, phe, withcov, regon=None, boost_modeld=None):
    # methodb='boosting_integrate'
    assert method == 'boosting'

    allstatd = {}
    for method_genet in ['boosting_add', 'boosting_nonadd']:
        # print(dresult, method, phe, withcov, regon, boost_modeld[method_genet])
        allstat = load_allstat(dresult, method, phe, withcov, regon, boost_modeld[method_genet])
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


def load_best_para(allstat, cvi, method, max_on, va_on, fcand,
                   return_stat=False):

    nsnvs, lrs = load_cand(fcand)

    # print('alltat in bestpara',allstat)

    allstat = allstat.loc[allstat['dataset'] == va_on, :]
    allstat = allstat.loc[allstat['cvi'] == cvi, :]

    if len(allstat) == 0:
        logger.warning('WARNING: dataset is not in allstat: {}, {}, {}'.format(method, cvi, va_on))
        return None

    # limit nsnvs
    if method == 'boosting':
        if nsnvs is not None:
            allstat = allstat.loc[allstat['n'].isin(nsnvs), :]
    elif method in iof.methods_classic_use_nsnvs():
        if nsnvs is not None:
            # allow all of more than largest
            max_nsnv = max(nsnvs)
            allstat = allstat.loc[allstat['n'].isin(nsnvs) | (allstat['n'] >= max_nsnv) | (allstat['n'] == -1), :]

    if method == 'boosting':
        assert 'lr' in allstat.columns
        # None : use the same as in config
        if (lrs != ()) and (lrs is not None):
            allstat = allstat.loc[allstat['lr'].isin(lrs), :]
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


def load_best_para_integrate(allstatd, cvi, method, max_on, va_on, fcand, return_stat=False):

    paras_best_add, stat_best_add = load_best_para(allstatd['boosting_add'], cvi, method, max_on, va_on, fcand, return_stat=True)
    # print('paras_best_add',paras_best_add)
    paras_best_nonadd, stat_best_nonadd = load_best_para(allstatd['boosting_nonadd'], cvi, method, max_on, va_on, fcand, return_stat=True)

    if stat_best_add > stat_best_nonadd:
        paras_best_add |= {'method_integrate': 'boosting_add'}
        if return_stat:
            return paras_best_add, stat_best_add
        else:
            return paras_best_add
    else:
        paras_best_nonadd |= {'method_integrate': 'boosting_nonadd'}
        if return_stat:
            return paras_best_nonadd, stat_best_nonadd
        else:
            return paras_best_nonadd


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

    # print('allstat', allstat)
    stat = allstat.loc[idx, stat_on]

    return stat


def load_best_stat_from_allstat(allstat, paras_best, dataset, stat_on, method, methodb, phe, max_on):
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
        'method': methodb,
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


def load_best_stat_from_allstat_integrate(allstatd, paras_best, dataset, stat_on, method, methodb, phe, max_on):
    if paras_best is None:
        return None

    method_integrate = paras_best['method_integrate']
    allstat = allstatd[method_integrate]

    statd = load_best_stat_from_allstat(allstat, paras_best, dataset, stat_on, method, methodb, phe, max_on)

    return statd


def load_best_stat_from_allstat_cvs_datasets(dresult, methodb, phe, withcov, reg_on, stat_on, max_on, cvn, boost_modeld, paras_best_cvs, datasets):
    # 'boosting_nonadd' -> 'boosting'
    method = methodb_to_method(methodb)

    stats = []
    if methodb == 'boosting_integrate':
        # allstat for stat
        allstatd = load_allstatd_integrate(dresult, method, phe, withcov, reg_on, boost_modeld)
        for cvi in range(cvn):
            paras_best = paras_best_cvs[cvi]
            for dataset in datasets:
                statd = load_best_stat_from_allstat_integrate(allstatd, paras_best, dataset, stat_on, method, methodb, phe, max_on)
                stats.append(statd)
    else:
        # None if not boosting
        boost_model = boost_modeld[methodb]

        # allstat for stat
        allstat = load_allstat(dresult, method, phe, withcov, reg_on, boost_model=boost_model)
        for cvi in range(cvn):
            paras_best = paras_best_cvs[cvi]
            for dataset in datasets:
                statd = load_best_stat_from_allstat(allstat, paras_best, dataset, stat_on, method, methodb, phe, max_on)
                stats.append(statd)
    return stats


def methodb_to_method(methodb):
    # if 'boosting_' in methodb:
    if methodb.startswith('boosting_'):
        method = 'boosting'
        # boost_model = boost_modeld[methodb]
    else:
        method = methodb
    return method


def load_best_para_cvs(dresult, methodb, phe, withcov_at_maxon, reg_on, max_on, va_on, fcand, boost_modeld, cvn=None, cvis=None):
    # 'boosting_nonadd' -> 'boosting'
    method = methodb_to_method(methodb)

    # None if not boosting
    boost_model = boost_modeld[methodb]

    # if 'boosting_' in methodb:
    #    method = 'boosting'
    #    boost_model = boost_modeld[methodb]
    # else:
    #    method = methodb
    #    boost_model = None

    # if methodb == 'boosting_integrate':
    #    boost_integrate = True
    # else:
    #    boost_integrate = False

    if cvn is not None:
        cvis = list(range(cvn))

    paras_best_cvs = {}
    if methodb == 'boosting_integrate':
        # allstat for maxon
        allstatd_maxon = load_allstatd_integrate(dresult, method, phe, withcov_at_maxon, reg_on, boost_modeld)
        # for cvi in range(cvn):
        for cvi in cvis:
            paras_best = load_best_para_integrate(allstatd_maxon, cvi, method, max_on, va_on, fcand)
            # print('paras_best',paras_best)
            # paras_best = load_best_para_integrate(allstatd_maxon, cvi, method, max_on, va_on, fymld['cand'])
            paras_best_cvs[cvi] = paras_best
    else:
        # allstat for maxon
        allstat_maxon = load_allstat(dresult, method, phe, withcov_at_maxon, reg_on, boost_model=boost_model)
        # for cvi in range(cvn):
        for cvi in cvis:
            paras_best = load_best_para(allstat_maxon, cvi, method, max_on, va_on, fcand)
            # paras_best = load_best_para(allstat_maxon, cvi, method, max_on, va_on, fymld['cand'])
            paras_best_cvs[cvi] = paras_best
    return paras_best_cvs


def regon_vaon(methodb, regond, withcov_at_maxon):
    # 'boosting_nonadd' -> 'boosting'
    method = methodb_to_method(methodb)

    if withcov_at_maxon:
        reg_on = regond[method]['regon']
    else:
        reg_on = None
    _, va_on = iof.get_datasets_vaon(reg_on)
    return reg_on, va_on


def get_method_model_para(methodb, boost_modeld, paras_best):
    if methodb == 'boosting_integrate':
        method_integrate = paras_best['method_integrate']
        boost_model = boost_modeld[method_integrate]
    else:
        # None if not boosting
        boost_model = boost_modeld[methodb]

    # 'boosting_nonadd' -> 'boosting'
    method = methodb_to_method(methodb)
    return (method, boost_model)


def extract_nsnv_best_para(paras_best_cvs, cvn, methodb, phe, max_on):
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
                'method': methodb,
                # 'method': method,
                'phe': phe,
                # 'dataset': dataset,
                # 'stat': stat,
                # 'stat_on': stat_on,
                'max_on': max_on,
                'cvi': paras_best['cvi'],
                # 'n': paras_best['n'],
                'nsnv_use': paras_best['nsnv_use'],
            }
        stats.append(statd)
    return stats


def load_wgt_str(fwgt):
    wgt = pd.read_csv(fwgt, delim_whitespace=True, dtype=str)
    return wgt


def extract_wgt_nsnv(wgt, nsnv):
    wgt.loc[:, 'nuniq'] = 0

    for rowi in range(len(wgt)):
        filt_snv = wgt['kind'] == 'SNV'
        nuniq = wgt.iloc[:rowi + 1, :].loc[filt_snv, 'var'].nunique()
        wgt.loc[rowi, 'nuniq'] = nuniq

    #print(wgt[['var', 'nuniq']].head(100).to_string())
    # row_nsnv = wgt[wgt['nuniq'] == nsnv].iloc[-1].index
    #print('row_nsnv', wgt[wgt['nuniq'] == nsnv])
    row_nsnv = wgt[wgt['nuniq'] == nsnv].index[-1]
    # row_nsnv = wgt[wgt['nuniq'] == nsnv].index.iloc[-1]
    #print('row_nsnv', row_nsnv)
    #print(wgt[['var', 'nuniq']].tail(100).to_string())

    # do not need +1 since .loc not .iloc
    wgt = wgt.loc[:row_nsnv, :]

    #print('wgt', wgt)

    return wgt


def phed_filename():
    phed = {
        't2d': 'type_2_diabetes',
        'cad': 'coronary_artery_disease',
        'bc': 'breast_cancer',
        'af': 'atrial_fibrillation',
        'ibd': 'inflammatory_bowel_disease',
        'gout': 'gout',
        'ra': 'rheumatoid_arthritis',
        'acd': 'all-cause_dementia',
        'ad': 'alzheimer\'s_disease',
        'psr': 'psoriasis',
        'cc': 'colorectal_cancer',
        'atm': 'asthma',
        'her': 'hernia',
        'ingher': 'inguinal_hernia',
        'ami': 'acute_myocardial_infarction',
        'mch': 'major_coronary_heart_disease_event',
        'mi': 'myocardial_infarction',
        'hae': 'hemorrhoids',
        'mbe': 'myopia_both_eyes'
    }
    return phed


def create_best_wgt(dout, dresult, methodbs, boost_modeld, phes, cvn, fymld, max_on, withcov_at_maxon):
    assert methodbs == ['boosting_integrate']

    regond = load_regon(fymld)
    # for methodb in methodbs:
    methodb = methodbs[0]
    reg_on, va_on = regon_vaon(methodb, regond, withcov_at_maxon)
    method = methodb_to_method(methodb)

    for phe in phes:
        paras_best_cvs = paper.load_best_para_cvs(dresult, methodb, phe, withcov_at_maxon, reg_on,
                                                  max_on, va_on, fymld['cand'], boost_modeld, cvn, )
        # nsnv_best_cvs = extract_nsnv_best_para(paras_best_cvs, cvn, methodb, phe, max_on)
        # stats.extend(nsnv_best_cvs)

        for cvi in range(cvn):
            paras_best = paras_best_cvs[cvi]

            method_integrate = paras_best['method_integrate']
            boost_model = boost_modeld[method_integrate]

            fwgt = iof.file_wgt(dresult, method, phe, cvi, paras_best, boost_model)
            wgt = load_wgt_str(fwgt)
            nsnv = paras_best['nsnv_use']
            #print('nsnv', nsnv)

            wgt = extract_wgt_nsnv(wgt, nsnv)
            #print('wgt', wgt.head())

            if 'score0' in wgt:
                cols = 'iteration,kind,var,model,alpha,const,score0,score1,score2,chrom,pos,a1,a2'.split(',')
            else:
                cols = 'iteration,kind,var,model,alpha,const,chrom,pos,a1,a2'.split(',')
            wgt = wgt[cols]

            phe_name = phed_filename()[phe]
            dphe = dout / phe_name
            os.makedirs(dphe, exist_ok=True)

            fout = dphe / ('cv' + str(cvi) + '.wgt')
            wgt.to_csv(fout, sep='\t', index=False, na_rep='NaN')

            if method_integrate == 'boosting_add':
                method_para = 'additive'
            elif method_integrate == 'boosting_nonadd':
                method_para = 'non-additive'
            else:
                raise RuntimeError
            dfpara = {
                'model': method_para,
                'learning-rate': paras_best['lr'],
                'num-snv': paras_best['nsnv_use']
            }
            s = pd.Series(dfpara)

            fout = dphe / ('cv' + str(cvi) + '.para')
            s.to_csv(fout, sep='\t', header=False)


def get_withcovs(max_on, plot_on=None,
                 # withcov_only=False,nocov_only=False, #ny
                 ):
    # TODO: not the best -> return list
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
        return [(True, True)]

    if max_on == 'nagelkerke':
        return [(False, True), (True, True)]

    if plot_on == 'nagelkerke':
        return [(True, False), (True, True)]

    return [(False, False), (False, True), (True, False), (True, True)]


def get_withcov_at_maxon(max_on):
    # withcov_at_maxon only
    # for nsnv

    if max_on == 'nagelkerke':
        return [True]

    return [False, True]


def dname_maxon(droot, max_on, withcov=None, withcov_at_maxon=None, phe=None, cvi=None):
    d = 'maxon_' + max_on

    if withcov_at_maxon is None:
        pass
    elif withcov_at_maxon:
        d += '.maxon_withcov'
    else:
        d += '.maxon_nocov'

    if withcov is None:
        pass
    elif withcov:
        d += '.withcov'
    else:
        d += '.nocov'

    if phe is None and cvi is None:
        d = droot / d
    else:
        if phe is not None and cvi is None:
            dir = phe
        elif phe is None and cvi is not None:
            dir = 'cv' + str(cvi)
        else:
            dir = phe + '.cv' + str(cvi)
        d = droot / d / dir

    # good idea
    # ds=['maxon_'+max_on,'maxonnocov_'+str(nocov_at_maxon),'nocov_'+str(nocov)]
    # dname='.'.join(ds)
    # d=droot+dname+'/'
    return d


def dname_maxon_create(droot, max_on, withcov=None, withcov_at_maxon=None, phe=None, cvi=None):
    d = dname_maxon(droot, max_on, withcov, withcov_at_maxon, phe, cvi)
    os.makedirs(d, exist_ok=True)
    return d


def create_public_wgt(dout, dresult, methodbs, boost_modeld, phes, cvn, fymld):

    for max_on in ['nagelkerke']:
        withcov_at_maxon = True
        create_best_wgt(dout, dresult, methodbs, boost_modeld, phes, cvn, fymld, max_on, withcov_at_maxon)


def extract_boosing_methods(methods, ext):
    methods_ext = []
    for method in methods:
        if method.startswith('boosting_') and method not in ext:
            continue
        else:
            methods_ext.append(method)
    return methods_ext


def public_wgt(dout, dresult, ddata, methodbs, boost_modeld,
               phes, cvn,
               para_cand, para_regon, para_com,
               ):

    fymld = get_fymld(ddata, para_cand, para_regon, para_com)

    create_public_wgt(dout, dresult, methodbs, boost_modeld, phes, cvn, fymld)


def main():
    args = argument()
    # 'boosting_integrate' not 'boosting'
    methodbs = args.methods.split(' ')
    phes = args.phes.split(' ')
    logger.debug('phes: {}'.format(phes))

    # 'boosting_nonadd' and 'boosting_add' are necessary for 'boostin

    boost_modeld = defaultdict(lambda: None)
    boost_modeld |= dict(args.method_model)
    # boost_modeld |= {
    #    'boosting_nonadd': args.model_boost,
    #    'boosting_add': args.model_boost_add,
    #    'boosting_loss': args.model_boost_loss
    # }
    logger.debug('boost_method {}'.format(boost_modeld))

    public_wgt(pathlib.Path(args.dout), pathlib.Path(args.dresult),
               pathlib.Path(args.ddata),
               methodbs, boost_modeld,
               phes, int(args.cvn),
               args.para_cand, args.para_regon, args.para_com,
               )


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
