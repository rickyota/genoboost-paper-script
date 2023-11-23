
import argparse
import pathlib
import numpy as np
from sklearn import metrics
import pandas as pd
import statsmodels.api as sm
# from collections import defaultdict
from logging import getLogger

from ...system.logger import logger_setting
from ...genetics.sample import sampleio as sampleio
# from ...genetics.io_genot import score as scoreio
from ...genetics.score import scoreio as scoreio
from ...genetics.sample import sample as sampleop
from ..io import filename as iof
from ..score import score as scoreop
from ..wgt import wgt as wgtop

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    # ./result/nonaddpgs/dataset.12/
    parser.add_argument('--dresult',
                        help='')
    #parser.add_argument('--dout',
                        #help='')

    parser.add_argument('--method',
                        help='')

    parser.add_argument('--mid-path',
                        nargs="+",
                        help='')

    parser.add_argument('--mid-path-covonly',
                        nargs="+",
                        help='')

    # parser.add_argument('--kind',
    #                    help='etc. default, nonadd, ex-chr6')

    # parser.add_argument('--model',
    #                    help='mainly for boosting')

    parser.add_argument('--regon',
                        help='')

    parser.add_argument('--regon_covonly',
                        help='')

    # parser.add_argument('--dscore',
    #                    help='')

    parser.add_argument('--cvn',
                        help='')

    parser.add_argument('--fcv',
                        help='')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--phe',
                        help='phenotype name')

    parser.add_argument('--sex',
                        help='')

    parser.add_argument('--mode',
                        help='')

    parser.add_argument('--nocov-only',
                        action='store_true',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


# def get_cv_names(cvs):
#    cols = cvs.columns
#    return [x for x in cols if x.startswith('cv')]


# def load_score(fin):
#    df = pd.read_csv(fin, header=0, delim_whitespace=True, dtype=str)
#    col_map_in = {'IID': 'id'}
#    df = df.rename(columns=col_map_in)
#    return df


# do not need para in allstat
# def create_allstat_boost(dscore, method, datasets, fcv):
#
#    fscore = dscore / 'boosting.score'
#    fpara = dscore / 'boosting.param'
#
#    score = load_score(fscore)
#
#    cvs = sampleio.load_fcv(fcv)
#    cv_names = get_cv_names(cvs)
#
#    for cv_name in cv_names:
#        for dataset in datasets:
#            #print("cv, dataset", cv_name, dataset)
#
#            sample_id = cvs.loc[cvs['cv0'] == dataset, 'id']
#
#            # untested
#            score_dataset = score.join(sample_id, on='id')
#
#            # acc
#
#    raise NotImplementedError()


# mv iof
# def paras_method(method):
#    if method == 'ldpred':
#        paras = ['rho', 'n']
#    elif method == 'clump':
#        paras = ['p', 'r2', 'n']
#    elif method == 'prscs':
#        paras = ['a', 'b', 'phi']
#    elif method == 'sbayesr':
#        paras = []
#    elif method == 'lassosum':
#        paras = ['s', 'n']
#    else:
#        raise NotImplementedError('Unknown method', method)
#    return paras


# mv io/sample.py
# def load_score_covonly(dout, phe, cvi, regon, mid_path, mode=None, reg_kind=None):
def load_score_covonly(dout, phe, cvi, regon, mid_path, mode=None, reg_kind=None):
    # mid_path = ('default',)
    fscorecov = iof.file_score(dout, 'covonly', phe, cvi, None, True, regon, mid_path, mode, reg_kind)
    logger.debug('fscore cov: {}'.format(fscorecov))
    score = scoreio.load_score(fscorecov)
    return score


# moved to scoreop
# def get_score_cols(score):
#    # or .startswith('score')
#    return [x for x in score.columns if (x != 'id') and (x != 'status')]
#    # return [x for x in score.columns if (x != 'id') and (x != 'phe')]


# moved to scoreop
# def get_n_para(col):
#    if col == 'score':
#        # return None
#        # use -1 only for pd => could forget
#        return -1
#    if col.startswith('score_n-'):
#        return int(col.replace('score_n-', ''))
#    raise RuntimeError('Cannot parse n_para: ', col)


# def get_nsnv_use(n_para, parad, dout, method, phe, cvi, kind, model):
def get_nsnv_use(n_para, parad, dout, method, phe, cvi, mid_path):
    if n_para != -1:
        return n_para

    # never in fscore
    # if 'n' in parad:
    #    return int(parad['n'])

    # load from wgt
    if method == 'covonly':
        return 0
    else:
        return wgtop.load_nsnv_wgt(dout, method, phe, cvi, parad, mid_path)
        # return wgtop.load_nsnv_wgt(dout, method, phe, cvi, parad, kind, model)


def calculate_pauc(pheno, fpr=None, score_col='score'):
    # fpr=None means usual auc
    # in order to distinguish not existing value, which is np.nan
    # pauc = metrics.roc_auc_score(phe['status'], phe['rs'], max_fpr=fpr)
    # pauc = metrics.roc_auc_score(pheno['status'], pheno[score_col], max_fpr=fpr)
    pauc = metrics.roc_auc_score(pheno['status'], pheno[score_col], max_fpr=fpr)

    # if np.isnan(pauc):
    #    pauc = 0.0

    return pauc


def calculate_auprc(pheno, score_col='score'):
    auprc = metrics.average_precision_score(pheno['status'], pheno[score_col])

    return auprc


def calc_likelihood_regression_score(rs, pheno):
    # llmodel= score + const
    # llnull = const
    const_in = np.ones_like(rs)
    X = np.stack([const_in, rs], axis=1)

    try:
        logit = sm.Logit(pheno, X).fit()

        # print(logit.summary())
        logger.debug('coef:\n{}'.format(logit.params))

        llf = logit.llf
        llnull = logit.llnull
    except np.linalg.LinAlgError as e:
        logger.info(e)
        llf, llnull = np.nan, np.nan
    except BaseException as e:
        logger.info('another error: ', e)
        # include 'PerfectSeparationError
        llf, llnull = np.nan, np.nan

    return llf, llnull


def calc_likelihood_regression(rs, pheno):
    # llmodel= genotonly + const
    # llnull = const

    return calc_likelihood_regression_score(rs, pheno)


def calc_likelihood_regression_withcov(rs, rs_cov, pheno):
    # llmodel= rs + const
    # llnull = rs_cov + const

    llmodel, _ = calc_likelihood_regression_score(rs, pheno)
    llnull, _ = calc_likelihood_regression_score(rs_cov, pheno)

    return llmodel, llnull


def calc_likelihood_metric(pheno, pheno_cov, withcov):
    # if nocov=F, model=phe and null=phe_cov
    # otherwise (T, 'va'), model=phe+const and null=const. (don't use phe_cov)

    if withcov:
        # TODO: ok when phe_cov include rs
        assert len(pheno['score']) == len(pheno['status']), 'length of rs and phe are different: ' + \
            str(len(pheno['score'])) + ' vs ' + str(len(pheno['status']))
        assert len(pheno['score']) == len(pheno_cov['score']), 'length of rs and rs cov are different: ' + \
            str(len(pheno['score'])) + ' vs ' + str(len(pheno_cov['score']))

        llmodel, llnull = calc_likelihood_regression_withcov(pheno['score'], pheno_cov['score'], pheno['status'])
        # loss funcs are different for {0,1} and {-1,1}
        # functions themselves are different but both output the exactly the same loss, so only use {0,1}
        # llmodel = calc_log_likelihood_model(phe['score'], phe['status'])
        # lmodel = np.exp(llmodel)
        logger.debug('llmodel: {}'.format(llmodel))

        # llnull = calc_log_likelihood_model(phe_cov['score'], phe['status'])
        # lnull = np.exp(llnull)
        logger.debug('llnull: {}'.format(llnull))
    else:
        assert len(pheno['score']) == len(pheno['status']), 'length of rs and phe are different: ' + \
            str(len(pheno['score'])) + ' vs ' + str(len(pheno['status']))

        llmodel, llnull = calc_likelihood_regression(pheno['score'], pheno['status'])

        # rs_reg = regression_nocov(phe['score'], phe['status'])
        # llmodel = calc_log_likelihood_model(rs_reg, phe['status'])
        logger.debug('llmodel: {}'.format(llmodel))
        # logger.debug('llmodel', llmodel)

        # rs_const = regression_constonly(phe['status'])
        # llnull = calc_log_likelihood_model(rs_const, phe['status'])
        logger.debug('llnull: {}'.format(llnull))
        # logger.debug('llnull', llnull)

    return llmodel, llnull


def calculate_nagelkerke(pheno, pheno_cov=None, withcov=True):

    assert set(pheno['status']) == {0, 1}

    # no 'status' in phe_cov
    # if withcov:
    #    assert set(phe_cov['status']) == {0, 1}

    llmodel, llnull = calc_likelihood_metric(pheno, pheno_cov, withcov)

    nsample = len(pheno)
    # nagel = (1 - (lnull / lmodel)**(2 / nsample)) / (1 - lnull**(2 / nsample))
    nagel = (1 - np.exp((llnull - llmodel) * (2 / nsample))) / (1 - np.exp(llnull * (2 / nsample)))
    # print('nagel', nagel)

    if nagel < 0.0 or 1.0 < nagel:
        logger.warning('**warning** nagel= {}'.format(nagel))
    # this might happen
    # assert 0.0 <= nagel <= 1.0

    return nagel


def calc_acc(score=None, score_covonly=None, withcov=None, baseline=False):
    accd = {}

    # auc
    if baseline:
        accd['auc'] = 0.5
    else:
        accd['auc'] = calculate_pauc(score)

    # pauc
    fprs = [0.03, 0.05, 0.08, 0.1, 0.3, 0.5]
    for fpr in fprs:
        acc_name = 'pauc_fpr-' + str(fpr)
        if baseline:
            accd[acc_name] = 0.5
        else:
            accd[acc_name] = calculate_pauc(score, fpr=fpr)

    # nagelkerke
    if baseline:
        accd['nagelkerke'] = 0.0
    else:
        accd['nagelkerke'] = calculate_nagelkerke(score, score_covonly, withcov)

    # auprc
    if baseline:
        # FIXME: depends on case prevalence
        accd['auprc'] = np.nan
    else:
        accd['auprc'] = calculate_auprc(score)

    return accd


def create_allstat_parad(dataset, cvi, parad, nsnv_use, accd):
    allstat_parad = {}
    allstat_parad |= {'dataset': dataset, 'cvi': cvi}
    allstat_parad |= parad
    allstat_parad |= {'nsnv_use': nsnv_use}
    allstat_parad |= accd
    # allstat_parad.update({'dataset': dataset, 'cvi': cvi})
    # allstat_parad.update(parad)
    # [allstat_parad.update({'n': n_para})
    # allstat_parad.update({'nsnv_use': nsnv_use})
    # allstat_parad.update(accd)
    return allstat_parad


def create_parad_col(parad, paras):
    # add '-1' if para not in parad
    parad_col = parad.copy()
    for para in paras:
        if para not in parad_col:
            parad_col[para] = '-1'
    return parad_col


def create_allstat(dresult, method, fcv, cvn, fpheno, phe, sex, withcov, regon, regon_covonly, mid_path, mid_path_covonly, mode):
    cvs = sampleio.load_cv(fcv)

    if 'dataset.12' in str(dresult) and method == 'snpnet':
        # TMP
        datasets = ['va', 'ts']
    else:
        # TODO: load from cvs
        datasets = ['tr', 'va', 'ts']

    allstat = []

    paras = iof.get_paras_method(method)

    if method == 'covonly' and (not withcov):
        # random prediction
        for cvi in range(cvn):
            for dataset in datasets:
                # if len(score_dataset) == 0:
                # continue
                accd = calc_acc(baseline=True)
                # accd = calc_acc(None, None, False, baseline=True)
                allstat_parad = create_allstat_parad(dataset, cvi, {}, 0, accd)
                # allstat_parad = create_allstat_parad(dataset, cvi, {}, -1, 0, accd)
                allstat.append(allstat_parad)
    else:
        for cvi in range(cvn):
            fscores = iof.files_score(dresult, method, phe, cvi, withcov, regon, mid_path, mode)
            # fscores = iof.fnames_score(dout, method, phe, cvi, withcov, regon, kind, model)
            logger.debug('fscores {}'.format(fscores))

            if withcov:
                # use score_covonly for method='covonly'??
                # -> yes, since no occasion to directly use score with cov
                # for auc, do not use score_covonly
                # for nagelkerke, covonly nagelkerke is not necessary

                score_covonly = load_score_covonly(dresult, phe, cvi, regon_covonly, mid_path_covonly, mode, reg_kind='logistic')
                logger.debug('score_covonly: {}'.format(score_covonly.head()))
                # score_covonly = load_score_covonly(dout, phe, cvi, 'tr', mode, reg_kind='logistic')
                # score_covonly = load_score_covonly(dout, phe, cvi, regon='tr', mid_path, mode, reg_kind='logistic')
            else:
                score_covonly = None

            for fscore in fscores:
                logger.debug('fscore {}'.format(fscore))
                parad = scoreop.parse_para(fscore)
                logger.debug('parad {}'.format(parad))
                score = scoreio.load_score_fpheno(fscore, fpheno, phe)
                cols = scoreop.get_score_cols(score)
                # print('cols', cols)
                # col='score' or 'score_n-100'
                for col in cols:
                    logger.debug('col {}'.format(col))
                    parad_col = create_parad_col(parad, paras)
                    # parad_col = parad.copy()

                    # get from col or wgt
                    # n as para
                    n_para = scoreop.get_n_para(col)
                    if 'n' in paras:
                        # even if n_para is -1, put as str
                        parad_col['n'] = str(n_para)
                    nsnv_use = get_nsnv_use(n_para, parad_col, dresult, method, phe, cvi, mid_path)
                    # nsnv_use = get_nsnv_use(n_para, parad_col, dout, method, phe, cvi, kind, model)
                    score_col = score.loc[:, ['id', 'status', col]].rename(columns={col: 'score'}).copy()
                    # nsnv_fname
                    for dataset in datasets:
                        score_dataset = sampleop.extract_dataset(score_col, dataset, cvi, cvs, sex)
                        if len(score_dataset) == 0:
                            continue
                        if withcov:
                            score_covonly_dataset = sampleop.extract_align(score_covonly, score_dataset)
                        else:
                            score_covonly_dataset = None
                        accd = calc_acc(score_dataset, score_covonly_dataset, withcov)

                        # update nsnv
                        allstat_parad = create_allstat_parad(dataset, cvi, parad_col, nsnv_use, accd)
                        allstat.append(allstat_parad)

    if len(allstat) == 0:
        cols = ['dataset', 'cvi'] + paras
        # cols = ['dataset', 'cvi'] + paras + ['n']
        # header only
        allstat = pd.DataFrame(columns=cols)
    else:
        allstat = pd.DataFrame.from_dict(allstat)
        # change column order
        cols_order = ['dataset', 'cvi'] + paras + ['nsnv_use']
        cols = cols_order + [x for x in allstat.columns.to_list() if x not in cols_order]
        allstat = allstat.loc[:, cols]

        sort_cols = ['dataset', 'cvi'] + paras
        allstat = allstat.sort_values(sort_cols).reset_index(drop=True)
    logger.debug('allstat {}'.format(allstat.head()))

    # TODO: pd.stack() ?
    return allstat


def create_allstats(dresult, method, mid_path, mid_path_covonly,
                    regon, regon_covonly,
                    fcv, cvn,
                    fpheno, phe, sex, mode, nocov_only):
    if nocov_only:
        withcovs = [False]
    else:
        withcovs = [True, False]

    for withcov in withcovs:

        allstat = create_allstat(dresult, method, fcv, cvn, fpheno, phe, sex, withcov, regon, regon_covonly, mid_path, mid_path_covonly, mode)

        fallstat = iof.file_allstat(dresult, method, phe, withcov, regon, mid_path, mode)
        logger.debug('fallstat {}'.format(fallstat))

        # no index
        allstat.to_csv(fallstat, sep='\t')
        # allstat.to_csv(fallstat, sep='\t',index=False)


def main():
    args = argument()

    # mid_path = [args.kind, args.model]

    if args.mid_path_covonly is None:
        mid_path_covonly = ('default',)
    else:
        mid_path_covonly = tuple(args.mid_path_covonly)

    if args.regon_covonly is None:
        regon_covonly = 'tr'
    elif args.regon_covonly in ['tr', 'va', 'ts']:
        regon_covonly = args.regon_covonly
    else:
        raise RuntimeError('Unknown regon_covonly: {}'.format(args.regon_covonly))

    create_allstats(pathlib.Path(args.dresult),
                    args.method,
                    # mid_path,
                    tuple(args.mid_path),
                    mid_path_covonly,
                    # tuple(args.mid_path_covonly),
                    # args.kind, args.model,
                    args.regon, regon_covonly,
                    args.fcv, int(args.cvn),
                    args.fpheno, args.phe,
                    args.sex, args.mode,
                    args.nocov_only)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
