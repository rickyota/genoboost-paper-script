
import os
import pathlib
import glob
from logging import getLogger


logger = getLogger(__name__)


def get_paras_method(method):
    if method == 'covonly':
        return []
    elif method == 'boosting':
        return ['lr', 'n']
    elif method == 'snpnet':
        return ['n']
    elif method == 'snpboost':
        # this works for no n score(i.e. snpboost.wgt)
        return ['n']
        # return []
    elif method == 'ldpred':
        return ['rho']
    elif method == 'clump':
        return ['p', 'r2', 'n']
    elif method == 'prscs':
        return ['a', 'b', 'phi']
    elif method == 'sbayesr':
        return []
    elif method == 'lassosum':
        return ['s', 'n']
    else:
        raise NotImplementedError('Unknown method', method)


def methods_classic_use_nsnvs():
    # TODO: add boosting?
    # for cand in allstat.py
    return ['clump', 'snpboost']


def _concat_paths(d, mid_path):
    if mid_path is not None:
        if not (isinstance(mid_path, list) or isinstance(mid_path, tuple)):
            raise RuntimeError('Use list or tuple for mid_path.')
        for d_path in mid_path:
            d /= d_path
    return d


# TODO: make phe in mid_path; not necessary since in paper.py phe is treated diferently
# def dir_result(dout, method, phe, kind=None, model=None, pattern=None):
def _dir_result(dresult, method, phe, mid_path=None):
    # dresult / mid_path[0] / mid_path[1] / ... / phe
    #
    # dresult / method / kind / pattern / model
    # TODO: mid_path=[kind, pattern, model] or mid=[kind, model] etc.
    # pattern is for simulation

    d = dresult / method

    d = _concat_paths(d, mid_path)

    d /= phe

    # TMP
    # if kind is None:
    #    raise NotImplementedError('Use kind="default".')

    # d = dout / method
    # if kind is not None:
    #    d /= kind
    # if pattern is not None:
    #    # Note pattern is middle
    #    d /= pattern
    # if model is not None:
    #    d /= model
    # d /= phe

    return d


def _is_in_parad(para, parad):
    # FIXME: what if true para could be -1
    # return ('n' in parad) and (parad['n'] != '-1' and parad['n'] != -1)
    return (para in parad) and (parad[para] != '-1' and parad[para] != -1)


# FIXME: make all parad string not int
# def fname_wgt_root(dout, method, phe, cvi, parad, kind=None, model=None, mode=None):
def _file_wgt_root(dresult, method, phe, cvi, parad, mid_path=None, mode=None, use_boosting_n=False):
    d = _dir_result(dresult, method, phe, mid_path)
    # d = dir_result(dout, method, phe, kind, model)

    if mode is None:
        # d = d / ('train.cv' + str(cvi))
        d /= ('train.cv' + str(cvi))
    else:
        # d = d / ('train_' + mode + '.cv' + str(cvi))
        d /= ('train_' + mode + '.cv' + str(cvi))

    if method == 'covonly':
        f = method
    elif method == 'boosting':
        f = method
        if _is_in_parad('lr', parad):
            f += '_lr-' + str(parad['lr'])
        if use_boosting_n and _is_in_parad('n', parad):
            # for regcov
            f += '_n-' + str(parad['n'])
    elif method == 'snpnet':
        f = method + '_n-' + str(parad['n'])
    elif method == 'snpboost':
        if _is_in_parad('n', parad):
            f = method + '_n-' + str(parad['n'])
        else:
            f = method
    elif method == 'ldpred':
        f = method + '_rho-' + str(parad['rho'])
    elif method == 'clump':
        # FIXME->ok?
        # if ('n' in parad) and (parad['n'] != -1):
        # hopefully parad should be string
        # if ('n' in parad) and (parad['n'] != '-1' and parad['n'] != -1):
        if _is_in_parad('n', parad):
            # if ('n' in parad) and (parad['n'] is not None) and (parad['n'] !=-1):
            f = method + '_p-' + str(parad['p']) + '_r2-' + str(parad['r2']) + '_n-' + str(parad['n'])
        else:
            f = method + '_p-' + str(parad['p']) + '_r2-' + str(parad['r2'])
    elif method == 'prscs':
        f = method + '_a-' + str(parad['a']) + '_b-' + str(parad['b']) + '_phi-' + str(parad['phi'])
    elif method == 'sbayesr':
        f = method
    elif method == 'lassosum':
        f = method + '_s-' + str(parad['s']) + '_n-' + str(parad['n'])
    else:
        raise NotImplementedError('Unknown method', method)
    return d, f


# def filename_wgt(dout, method, phe, cvi, parad, kind=None, model=None, mode=None):
def file_wgt(dresult, method, phe, cvi, parad, mid_path=None, mode=None):
    # parad should have string value
    d, f = _file_wgt_root(dresult, method, phe, cvi, parad, mid_path, mode)
    # d, f = fname_wgt_root(dout, method, phe, cvi, parad, kind, model, mode)
    fwgt = d / (f + '.wgt')
    return fwgt


def _file_wgt_boosting_root_best_para(dresult, method, phe, cvi, mid_path, mode=None):
    d = _dir_result(dresult, method, phe, mid_path)

    if mode is None:
        d /= ('train.cv' + str(cvi))
    else:
        d /= ('train_' + mode + '.cv' + str(cvi))

    if method == 'boosting':
        f = method
    else:
        raise NotImplementedError('Unknown method', method)
    return d, f


def file_wgt_boosting_best_para(dresult, method, phe, cvi, mid_path, mode=None):
    d, f = _file_wgt_boosting_root_best_para(dresult, method, phe, cvi, mid_path, mode=mode)
    fwgt = d / (f + '.wgt')
    return fwgt


# merge above: how to flexibly add layer?
# def filename_wgt_boosting_best_para_simu(dout, method, phe, cvi, mid_path, mode=None):
#    print(dout, method, phe, cvi, kind, model, simu_type, mode)
#    d, f = fname_wgt_boosting_root_best_para(dout, method, phe, cvi, kind, model, simu_type, mode)
#    fwgt = d / (f + '.wgt')
#    return fwgt


# parad should have string value
def file_loss(dresult, phe, cvi, lr, iteri, mid_path):
    method = 'boosting'
    d = _dir_result(dresult, method, phe, mid_path)
    d = d / ('train.cv' + str(cvi)) / 'loss' / ('para_lr-' + lr)
    # error since phe is str
    # d /= phe / ('train.cv' + str(cvi)) / 'loss' / ('para_lr-' + lr)
    # f = d / phe / ('train.cv' + str(cvi)) / 'loss' / ('para_lr-' + lr)
    # d = dout / method / kind / model / phe / ('train.cv' + str(cvi)) / 'loss' / ('para_lr-' + lr)
    f = 'iter-' + str(iteri) + '.loss'
    return d / f


# merge to fname_wgtcov()
# def fname_wgtcov_covonly(dout, phe, cvi, regon, reg_kind):
#    if reg_kind == 'logistic':
#        return dout / 'covonly' / phe / ('train.cv' + str(cvi)) / ('covonly.wgtcov_regon' + regon + '.gz')
#    else:
#        raise NotImplementedError('reg_kind: ', reg_kind)
#    # return dout / ('covonly.wgtcov_regon' + regon)


def file_wgtcov(dresult, method, phe, cvi, parad, regon, mid_path, mode=None, reg_kind=None):
    # parad should have string value
    d, f = _file_wgt_root(dresult, method, phe, cvi, parad, mid_path, mode, use_boosting_n=True)
    # print('fwgt_root', d, f)
    if regon is None:
        ext = '.wgtcov'
    else:
        ext = '.wgtcov_regon' + regon

    # TODO: reg_kind -> model ?
    if (reg_kind is None) or (reg_kind == 'logistic'):
        # mainly for covonly
        pass
    else:
        raise NotImplementedError('reg_kind: ', reg_kind)

    fwgtcov = d / (f + ext)
    return fwgtcov


# TODO: reg_kind -> model ?
# def fname_scorecov_covonly(dresult, phe, cvi, regon, reg_kind):
#    if reg_kind == 'logistic':
#        return dresult / 'covonly' / 'default' / phe / ('score.cv' + str(cvi)) / ('covonly.scorecov_regon' + regon + '.gz')
#    else:
#        raise NotImplementedError('reg_kind: ', reg_kind)


def file_score(dresult, method, phe, cvi, parad, withcov=False, regon=None, mid_path=None, mode=None, reg_kind=None):

    d = _dir_result(dresult, method, phe, mid_path)

    if mode is None:
        d /= ('score.cv' + str(cvi))
    else:
        d /= ('score_' + mode + '.cv' + str(cvi))

    if method == 'covonly':
        f = method
    elif method == 'boosting':
        # if 'lr' in parad:
        if _is_in_parad('lr', parad):
            f = method + '_lr-' + str(parad['lr']) + '_n'
            # f = method + '_lr-' + parad['lr'] + '_n' + parad['n'] + '.score'
        else:
            f = method + '_n'
    elif method == 'snpnet':
        f = method + '_n'
    elif method == 'snpboost':
        if _is_in_parad('n', parad):
            f = method + '_n'
        else:
            f = method
    elif method == 'ldpred':
        f = method + '_rho-' + str(parad['rho'])
    elif method == 'clump':
        # FIXME:ok? use int?
        # if ('n' in parad) and (parad['n'] != -1):
        if _is_in_parad('n', parad):
            # if ('n' in parad) and (parad['n'] is not None) and (parad['n'] =-1):
            f = method + '_p-' + str(parad['p']) + '_r2-' + str(parad['r2']) + '_n'
        else:
            f = method + '_p-' + str(parad['p']) + '_r2-' + str(parad['r2'])
    elif method == 'prscs':
        f = method + '_a-' + str(parad['a']) + '_b-' + str(parad['b']) + '_phi-' + str(parad['phi'])
    elif method == 'sbayesr':
        f = method
    elif method == 'lassosum':
        f = method + '_s-' + str(parad['s']) + '_n'
    else:
        raise NotImplementedError('Unknown method', method)

    # withcov
    if method == 'covonly' and (not withcov):
        raise RuntimeError('Usually use withcov=True')

    if withcov:
        if regon is None:
            f += '.scorecov'
        else:
            f += '.scorecov_regon' + regon
    else:
        f += '.score'

    f += '.gz'

    # TODO: reg_kind -> model ?
    if (reg_kind is None) or (reg_kind == 'logistic'):
        # mainly for covonly
        pass
    else:
        raise NotImplementedError('reg_kind: ', reg_kind)

    return d / f


def files_score(dout, method, phe, cvi, withcov, regon, mid_path, mode=None):
    """
    mode: ex. 'split-add'
    """

    d = _dir_result(dout, method, phe, mid_path)

    if mode is None:
        d = d / ('score.cv' + str(cvi))
    else:
        d = d / ('score_' + mode + '.cv' + str(cvi))

    if withcov:
        if regon is None:
            # for snpnet
            f = method + '*.scorecov'
        else:
            f = method + '*.scorecov_regon' + regon
    else:
        f = method + '*.score'
    f += '.gz'

    g_fscore = d / f
    logger.debug('g_fscore {}'.format(g_fscore))
    fscores = glob.glob(str(g_fscore))
    fscores = [pathlib.Path(f) for f in fscores if os.path.isfile(f)]
    # sometime, when not exist, return '*.score' itself
    # fscores = [pathlib.Path(f) for f in fscores]
    return fscores


def file_allstat(dresult, method, phe, withcov, regon, mid_path, mode=None):
    d = _dir_result(dresult, method, phe, mid_path)

    f = method + '.allstat'

    # if withcov:
    #    if regon is None:
    #        f = (method + '.allstat_cov')
    #    else:
    #        f = (method + '.allstat_cov_regon' + regon)
    # else:
    #    f = (method + '.allstat')

    if mode is not None:
        f += '_' + mode

    if withcov:
        if regon is None:
            f += '_cov'
        else:
            f += '_cov_regon' + regon

    fallstat = d / f
    # fallstat = dout / method / phe / stem
    return fallstat


def get_datasets_vaon(reg_on):
    if reg_on is None:
        datasets = ['va', 'ts']
        vaon = 'va'
    elif reg_on == 'tr':
        datasets = ['va', 'ts']
        vaon = 'va'
    elif reg_on == 'va':
        datasets = ['va', 'ts']
        vaon = 'va'
    elif reg_on == 'va1':
        raise RuntimeError('use regon_ori=tr')
        datasets = ['va1', 'va2', 'test']
        vaon = 'va2'
    else:
        raise NotImplementedError('ny', reg_on)
    return datasets, vaon


# parad should have string value
def file_ss(dresult, phe, cvi, gmodel='add', mid_path=None):

    d = dresult / 'assoc'

    d = _concat_paths(d, mid_path)

    # if mid_path is not None:
    #    if not isinstance(mid_path, list):
    #        raise RuntimeError('Use list for mid_path.')
    #    for d_path in mid_path:
    #        d /= d_path

    if gmodel == 'add':
        f_gmodel = ''
    else:
        f_gmodel = '.' + gmodel

    f = 'ukb_imp.cv' + str(cvi) + '.' + phe + f_gmodel + '.glm.logistic.hybrid.ss'

    return d / f
