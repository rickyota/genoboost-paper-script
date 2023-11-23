
import os
import argparse
import pathlib
import pandas as pd
from sklearn import metrics
import numpy as np
import statsmodels.api as sm
from logging import getLogger

from ...system.logger import logger_setting
from ...genetics.sample import sampleio as sampleio
from ...genetics.score import scoreio as scoreio
from ...genetics.sample import sample as sampleop
from ..io import filename as iof
# from ..io import arg as ioarg
from ..score import score as scoreop

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    # ./result/nonaddpgs/dataset.12/
    parser.add_argument('--dresult',
                        help='result dir')
    # parser.add_argument('--dout',
    #                    help='output dir')

    parser.add_argument('--train',
                        action='store_true',
                        help='for regression')

    parser.add_argument('--score',
                        action='store_true',
                        help='for score')

    parser.add_argument('--method',
                        help='ex. covonly, ')

    # parser.add_argument('--kind',
    #                    help='nonadd')

    # parser.add_argument('--model',
    #                    help='mainly for boosting')

    parser.add_argument('--mid-path',
                        nargs="+",
                        help='')

    parser.add_argument('--mode-score',
                        help='')

    parser.add_argument('--mode-train',
                        help='')

    parser.add_argument('--regon',
                        # required=True,
                        # None for snpnet
                        default=None,
                        help='ex. tr, va')

    parser.add_argument('--fcv',
                        help='')

    parser.add_argument('--cvi',
                        help='')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--phe',
                        help='phe')

    # parser.add_argument('--covs',
    parser.add_argument('--cov-names',
                        # nargs='+',
                        help='ex. age,sex,PC_1-PC_10. When use "-", the cov ("PC") should be the same.')

    parser.add_argument('--sex',
                        help='phe')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def regression_pheno(pheno, cov):
    # phe=phe_va.set_index('id')
    # cov=cov.set_index('id')

    # X_pd = pd.concat([phe_va['status'], cov], axis=1, join='inner')
    # logger.debug("pheno\n{}".format(pheno))
    # logger.debug("cov\n{}".format(cov))
    x_pd = pd.merge(pheno, cov, how='left', on='id')
    # logger.debug("x_pd\n{}".format(x_pd))

    x = x_pd.loc[:, [(x not in ['id', 'status']) for x in x_pd.columns]]
    # logger.debug("x\n{}".format(x.info()))
    # logger.debug("x\n{}".format(x))
    # logger.debug("x\n{}".format(x.columns))
    y = x_pd['status']
    # print('sety', set(y))

    assert set(y) == {0, 1}, 'otherwise sm.Logit would output error'

    # using statsmodels for 95% confidential interval
    x = sm.add_constant(x)
    # print("add const", x.shape)

    logit = sm.Logit(y, x).fit()

    logger.debug('{}'.format(logit.summary()))

    logit_summary = pd.read_html(logit.summary().tables[1].as_html(), header=0, index_col=0)[0]
    return logit_summary


# def fname_wgtcov_covonly(dout, phe_name, cvi, regon, reg_kind):
#    if reg_kind == 'logistic':
#        return dout / phe_name / ('tr.cv' + str(cvi)) / ('covonly.wgtcov_regon' + regon)
#    else:
#        raise NotImplementedError('reg_kind: ', reg_kind)
#    # return dout / ('covonly.wgtcov_regon' + regon)


# def fname_scorecov_covonly(dout, phe_name, cvi, regon, reg_kind):
#    if reg_kind == 'logistic':
#        return dout / phe_name / ('score.cv' + str(cvi)) / ('covonly.scorecov_regon' + regon)
#    else:
#        raise NotImplementedError('reg_kind: ', reg_kind)
#    # return dout / ('covonly.wgtcov_regon' + regon)


def regression_withcov_covonly(dresult,
                               regon,
                               fcv, cvi,
                               fpheno, phe, covs,
                               sex,
                               mid_path, mode_train,
                               reg_kind='logistic'):
    method = 'covonly'

    covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
    logger.debug('cov\n{}'.format(covar.head()))

    pheno = sampleio.load_pheno_phe(fpheno, phe)
    logger.debug('phe\n{}'.format(pheno.head()))
    cvs = sampleio.load_cv(fcv)
    logger.debug('cvs\n{}'.format(cvs.head()))

    pheno_ext = sampleop.extract_dataset(pheno, dataset=regon, cvi=cvi, cvs=cvs, sex=sex)
    logger.debug('phe_va\n{}'.format(pheno_ext.head()))
    if len(pheno_ext) == 0:
        raise RuntimeError('samples=0')
    # not necessary but ok
    covar_ext = sampleop.extract_align(covar, pheno_ext)
    logger.debug('cov_va\n{}'.format(covar_ext.head()))

    if reg_kind == "logistic":
        try:
            logit_summary = regression_pheno(pheno_ext, covar_ext)
        except np.linalg.LinAlgError as e:
            logger.error('LinAlgError: {}'.format(e))
            return

        fwgtcov = iof.file_wgtcov(dresult, method, phe, cvi, None, regon, mid_path, mode_train, reg_kind)
        # fwgtcov = iof.fname_wgtcov_covonly(dout, phe, cvi, regon, reg_kind)
        make_dir_file(fwgtcov)
        logit_summary.to_csv(fwgtcov, sep='\t', na_rep='NaN')
    else:
        raise RuntimeError('Unknown reg_kind:', reg_kind)


def regression_withcov_classic(dout,
                               method,
                               regon,
                               fcv, cvi,
                               fpheno, phe, covs,
                               sex, mid_path, mode_train, reg_kind='logistic'):

    covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
    logger.debug('cov\n{}'.format(covar.head()))

    pheno = sampleio.load_pheno_phe(fpheno, phe)
    logger.debug('phe\n{}'.format(pheno.head()))
    cvs = sampleio.load_cv(fcv)
    logger.debug('cvs\n{}'.format(cvs.head()))

    # pheno_ext = sampleop.extract_dataset(pheno, dataset=regon, cvi=cvi, cvs=cvs, sex=sex)
    # logger.debug('phe_va\n{}'.format(pheno_ext.head()))
    # not necessary but ok
    # covar_ext = sampleop.extract_align(covar, pheno_ext)
    # logger.debug('cov_va\n{}'.format(covar_ext.head()))

    if reg_kind == "logistic":

        fscores = iof.files_score(dout, method, phe, cvi, False, regon, mid_path)
        for fscore in fscores:
            parad = scoreop.parse_para(fscore)

            score_ori = scoreio.load_score_fpheno(fscore, fpheno, phe)
            score = sampleop.extract_dataset(score_ori, dataset=regon, cvi=cvi, cvs=cvs, sex=sex)
            if len(score) == 0:
                raise RuntimeError('#samples=0.')

            cols = scoreop.get_score_cols(score)

            for col in cols:
                # print('col', col)
                parad_col = parad.copy()
                n_para = scoreop.get_n_para(col)
                if n_para != -1:
                    # update parad
                    parad_col['n'] = str(n_para)
                score_col = score.loc[:, ['id', 'status', col]].rename(columns={col: 'score'}).copy()

                fwgtcov = iof.file_wgtcov(dout, method, phe, cvi, parad_col, regon, mid_path)
                # print('fwgtcov', fwgtcov)

                try:
                    logit_summary = regression_pheno(score_col, covar)
                except np.linalg.LinAlgError as e:
                    logger.error('LinAlgError: {}'.format(e))
                    # print(e)
                    return

                # print('logitsum', logit_summary)

                make_dir_file(fwgtcov)
                logit_summary.to_csv(fwgtcov, sep='\t', na_rep='NaN')
    else:
        raise RuntimeError('Unknown reg_kind:', reg_kind)


def addcov_score(pheno, covar, logitsum):
    x_pd = pd.merge(pheno, covar, how='outer', on='id')
    x_pd.loc[:, 'const'] = 1.0
    # logger.debug("x_pd\n{}".format(x_pd.info()))
    # logger.debug("x_pd\n{}".format(x_pd.columns))
    # logger.debug("x_pd\n{}".format(x_pd.head()))

    # for snpnet
    if 'score' not in logitsum.index:
        logitsum.loc['score', 'coef'] = 1.0

    cols = logitsum.index
    if set(x_pd.columns) < set(cols):
        raise RuntimeError("Some columns in logitsum are not in x: {}, {}", x_pd.columns, "vs", logitsum.index)
    # logger.debug("logitsum\n{}".format(logitsum.index))
    # logger.debug("logitsum\n{}".format(logitsum.head()))
    x = x_pd.loc[:, cols]
    # logger.debug("x\n{}".format(x.columns))
    # logger.debug("x\n{}".format(x.head()))

    if (x.columns != logitsum.index).any():
        raise RuntimeError("Cov names are different in x and logitsum: ", x.columns, "vs", logitsum.index)

    y_pred = x.values @ logitsum['coef'].values

    # just for log
    # auc = metrics.roc_auc_score(pheno['status'], y_pred)
    # logger.debug('auc direct {}'.format(auc))

    pheno_score = pheno[['id']].copy()
    pheno_score.loc[:, 'score'] = y_pred

    return pheno_score


# could be same as addcov_score?
def addcov_score_covonly(pheno, covar, logitsum):
    x_pd = pd.merge(pheno, covar, how='outer', on='id')
    x_pd.loc[:, 'const'] = 1.0
    cols = logitsum.index
    logger.debug("logitsum\n{}".format(logitsum.info()))
    x = x_pd.loc[:, cols]
    logger.debug("x\n{}".format(x.info()))

    # if (x.columns != logitsum.index).all():
    if (x.columns != logitsum.index).any():
        raise RuntimeError("Cov names are different in x and logitsum: ", x.columns, "vs", logitsum.index)

    y_pred = x.values @ logitsum['coef'].values

    # just for log
    auc = metrics.roc_auc_score(pheno['status'], y_pred)
    logger.debug('auc direct: {}'.format(auc))

    pheno_score = pheno[['id']].copy()
    pheno_score.loc[:, 'score'] = y_pred

    return pheno_score


def make_dir_file(f):
    os.makedirs(f.parent, exist_ok=True)
    # os.makedirs(pathlib.Path(f).parent)


def score_withcov_covonly(dresult,
                          regon,
                          fcv, cvi,
                          fpheno, phe, covs,
                          # sex,  # unnecessary
                          mid_path,
                          mode_score, mode_train,
                          reg_kind='logistic'):
    method = 'covonly'

    covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
    logger.debug('cov\n{}'.format(covar.head()))

    pheno = sampleio.load_pheno_phe(fpheno, phe)
    logger.debug('phe\n{}'.format(pheno.head()))
    cvs = sampleio.load_cv(fcv)
    logger.debug('cvs\n{}'.format(cvs.head()))

    # do not extract dataset; calculate all scores
    # phe_ext = sampleop.extract_dataset(phe, dataset=regon, cvs=cvs, cvi=cvi, sex=sex)
    # logger.debug('phe_va\n{}'.format(phe_ext.head()))
    # cov_va = sampleop.extract_align(cov, phe_ext)
    # logger.debug('cov_va\n{}'.format(cov_va.head()))

    if reg_kind == 'logistic':
        fwgtcov = iof.file_wgtcov(dresult, method, phe, cvi, None, regon, mid_path, mode_train, reg_kind)
        # fwgtcov = iof.fname_wgtcov_covonly(dresult, phe, cvi, regon, reg_kind)
        logit_summary = pd.read_csv(fwgtcov, sep='\t', index_col=0)

        score_covonly = addcov_score_covonly(pheno, covar, logit_summary)

        fscorecov = iof.file_score(dresult, method, phe, cvi, None, True, regon, mid_path, mode_score, reg_kind)
        make_dir_file(fscorecov)
        # score_covonly.to_csv(fscorecov, sep='\t', index=False)
        scoreio.write_score(fscorecov, score_covonly)

    else:
        raise RuntimeError('Unknown reg_kind:', reg_kind)


def score_withcov_classic(dresult,
                          method,
                          regon,
                          # fcv,
                          cvi,
                          fpheno, phe,
                          covs,
                          # kind, model,
                          mid_path,
                          mode_score, mode_train,
                          reg_kind='logistic'):

    covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
    # covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
    logger.debug('cov\n{}'.format(covar.head()))

    # pheno = sampleio.load_pheno_phe(fpheno, phe)
    # logger.debug('pheno\n{}'.format(pheno.head()))
    # calculate all samples
    # cvs = sampleio.load_cv(fcv)
    # logger.debug('cvs\n{}'.format(cvs.head()))

    # method = 'snpnet'

    if reg_kind == 'logistic':
        fscores = iof.files_score(dresult, method, phe, cvi, False, regon, mid_path, mode_score)

        # print('fscores', fscores)
        for fscore in fscores:
            logger.debug('fscore: {}'.format(fscore))
            parad = scoreop.parse_para(fscore)

            score = scoreio.load_score_fpheno(fscore, fpheno, phe)
            # print('score', score.head())

            cols = scoreop.get_score_cols(score)

            score_cols = []
            for col in cols:
                parad_col = parad.copy()
                n_para = scoreop.get_n_para(col)
                if n_para != -1:
                    # update parad
                    parad_col['n'] = str(n_para)
                score_col = score.loc[:, ['id', 'status', col]].rename(columns={col: 'score'}).copy()

                fwgtcov = iof.file_wgtcov(dresult, method, phe, cvi, parad_col, regon, mid_path, mode_train)
                logit_summary = pd.read_csv(fwgtcov, sep='\t', index_col=0)

                score_col_withcov = addcov_score(score_col, covar, logit_summary)
                score_cols.append(score_col_withcov)

            score_withcov = score_cols[0].loc[:, ['id']].copy()
            for df, col in zip(score_cols, cols):
                score_withcov = pd.merge(score_withcov, df.rename(columns={'score': col}), on='id', how='outer')
            logger.debug("score_withcov\n{}".format(score_withcov.head()))

            # FIXME: cannot distinguish '_n' or not
            # problem for clump and snpboost
            # TMP
            if str(fscore).endswith('_n.score.gz'):
                parad['n'] = '1'
            fscore_withcov = iof.file_score(dresult, method, phe, cvi, parad, True, regon, mid_path, mode_score)
            logger.debug('fscore_withcov: {}'.format(fscore_withcov))
            make_dir_file(fscore_withcov)
            # score_covonly.to_csv(fscorecov, sep='\t', index=False)
            scoreio.write_score(fscore_withcov, score_withcov)

    else:
        raise RuntimeError('Unknown reg_kind:', reg_kind)

# def score_withcov_snpnet(dout,
#                         regon,
#                         fcv, cvi,
#                         fpheno, phe,
#                         covs,
#                         model,
#                         reg_kind='logistic'):
#
#    covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
#    # covar = sampleio.load_cov(fpheno, cov_names_unparse=covs)
#    logger.debug('cov\n{}'.format(covar.head()))
#
#    pheno = sampleio.load_pheno_phe(fpheno, phe)
#    logger.debug('phe\n{}'.format(pheno.head()))
#    cvs = sampleio.load_cv(fcv)
#    logger.debug('cvs\n{}'.format(cvs.head()))
#
#    method = 'snpnet'
#
#    if reg_kind == 'logistic':
#        fscores = iof.fnames_score(dout, method, phe, cvi, False, regon, model)
#
#        # print('fscores', fscores)
#        for fscore in fscores:
#            parad = scoreop.parse_para(fscore)
#
#            score = scoreio.load_score_fpheno(fscore, fpheno, phe)
#            # print('score', score.head())
#
#            cols = scoreop.get_score_cols(score)
#
#            score_cols = []
#            for col in cols:
#                parad_col = parad.copy()
#                n_para = scoreop.get_n_para(col)
#                if n_para != -1:
#                    # update parad
#                    parad_col['n'] = str(n_para)
#                score_col = score.loc[:, ['id', 'status', col]].rename(columns={col: 'score'}).copy()
#
#                fwgtcov = iof.fname_wgtcov(dout, method, phe, cvi, parad_col, regon, model)
#                # fwgtcov = iof.fname_wgtcov_snpnet(dout, phe, cvi, parad_col, regon, model)
#                logit_summary = pd.read_csv(fwgtcov, sep='\t', index_col=0)
#
#                # fscore_va = iof.filename_score(dout, method, phe, cvi, parad, False, regon, model)
#                # score = scoreio.load_score_fpheno(fscore_va, fpheno, phe)
#
#                score_col_withcov = addcov_score(score_col, covar, logit_summary)
#                # print('score_col_withcov', score_col_withcov.head())
#                # raise NotImplementedError('ny')
#                score_cols.append(score_col_withcov)
#
#            # score_cols=[df.set_index('id') for df in score_cols]
#            score_withcov = score_cols[0].loc[:, ['id']].copy()
#            for df, col in zip(score_cols, cols):
#                score_withcov = pd.merge(score_withcov, df.rename(columns={'score': col}), on='id', how='outer')
#            logger.debug("score_withcov\n{}".format(score_withcov.head()))
#
#            fscore_withcov = iof.filename_score(dout, method, phe, cvi, parad, True, regon, model)
#            make_dir_file(fscore_withcov)
#            # score_covonly.to_csv(fscorecov, sep='\t', index=False)
#            scoreio.write_score(fscore_withcov, score_withcov)
#
#    else:
#        raise RuntimeError('Unknown reg_kind:', reg_kind)


def main():
    args = argument()

    # mid_path = args.mid_path
    # logger.debug('mid_path {}'.format(mid_path))

    method = args.method

    if method == 'covonly':
        if args.train:
            if args.mode_train is not None:
                raise RuntimeError('mode not implemented')
            regression_withcov_covonly(pathlib.Path(args.dresult),
                                       args.regon,
                                       args.fcv, args.cvi,
                                       args.fpheno, args.phe, args.cov_names,
                                       args.sex,
                                       tuple(args.mid_path), args.mode_train
                                       )
        if args.score:
            score_withcov_covonly(pathlib.Path(args.dresult),
                                  args.regon,
                                  args.fcv, args.cvi,
                                  args.fpheno, args.phe, args.cov_names,
                                  # args.sex,
                                  tuple(args.mid_path),
                                  args.mode_score, args.mode_train)
    # elif method == 'snpnet':
    #    if args.train:
    #        raise RuntimeError('No training here.')
    #    if args.score:
    #        score_withcov_snpnet(pathlib.Path(args.dresult),
    #                             args.regon,
    #                             args.fcv, args.cvi,
    #                             args.fpheno, args.phe, args.cov_names, args.model)
    else:
        if args.train:
            if args.regon == 'tr' and (method == 'snpnet' or method == 'snpboost'):
                raise RuntimeError('No training here.')
            regression_withcov_classic(pathlib.Path(args.dresult),
                                       args.method,
                                       args.regon,
                                       args.fcv, args.cvi,
                                       args.fpheno, args.phe, args.cov_names, args.sex,
                                       tuple(args.mid_path), args.mode_train)
            # args.kind, args.model)
        if args.score:
            score_withcov_classic(pathlib.Path(args.dresult),
                                  args.method,
                                  args.regon,
                                  # args.fcv,
                                  args.cvi,
                                  args.fpheno, args.phe, args.cov_names,
                                  tuple(args.mid_path), args.mode_score, args.mode_train)
            # args.kind, args.model)
    # else:
    #    raise RuntimeError('Unknown method: ', args.method)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
