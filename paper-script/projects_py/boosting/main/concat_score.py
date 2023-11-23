
"""
concatenate scores by aggreagating #snvs
"""

import argparse
from logging import getLogger
from ...system.logger import logger_setting
import pandas as pd
import pathlib
import glob


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--dout',
                        help='output file')

    parser.add_argument('--dout-score',
                        help='output file')

    parser.add_argument('--concat-paras',
                        action='store_true',
                        help='')

    parser.add_argument('--method',
                        help='')

    parser.add_argument('--para-concat',
                        help='')

    parser.add_argument('--regon',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


# this should be faster and more accurate than using regex
def parse_n(f):
    # assume suffix only has one '.'
    # ok: '.score_withcov'
    # ng: '.score.withcov'
    base = pathlib.Path(f).stem
    pnames = base.split('_')
    for pname in pnames:
        if '-' not in pname:
            continue
        ps = pname.split('-')
        if ps[0] != 'n':
            continue
        # para of n
        return ps[1]
    raise RuntimeError('Parameter n is not in filename: ', f)


def get_fscores_n(d, para_fix, para, method, regon):
    assert para == 'n'
    if method == 'boosting':
        # regon is same since d contains withcov
        g_fscore = d / (method + '_n-*.score')
    elif method == 'snpnet':
        if regon is None:
            g_fscore = d / (method + '_n-*.score')
        elif regon == 'tr':
            # TMP
            g_fscore = d / (method + '_n-*.scorecov')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'clump':
        if regon is None:
            g_fscore = d / ('clump_p-' + para_fix['p'] + '_r2-' + para_fix['r2'] + '_n-*.score')
        elif regon == 'va':
            g_fscore = d / ('clump_p-' + para_fix['p'] + '_r2-' + para_fix['r2'] + '_n-*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'lassosum':
        if regon is None:
            g_fscore = d / ('lassosum_s-' + para_fix['s'] + '_n-*.score')
        elif regon == 'va':
            g_fscore = d / ('lassosum_s-' + para_fix['s'] + '_n-*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    else:
        raise NotImplementedError('Unknown method', method)
    fscores = glob.glob(str(g_fscore))

    nsnvs = [parse_n(f) for f in fscores]

    return fscores, nsnvs


def get_fscore_concat(d, para_fix, para_concat, method, regon):
    assert para_concat == 'n'
    if method == 'boosting':
        fscore = d / (method + '_n.score')
    elif method == 'snpnet':
        if regon is None:
            fscore = d / (method + '_n.score')
        elif regon == 'tr':
            # TMP
            fscore = d / (method + '_n.scorecov_regontr')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'clump':
        if regon is None:
            fscore = d / ('clump_p-' + para_fix['p'] + '_r2-' + para_fix['r2'] + '_n.score')
        elif regon == 'va':
            fscore = d / ('clump_p-' + para_fix['p'] + '_r2-' + para_fix['r2'] + '_n.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'lassosum':
        if regon is None:
            fscore = d / ('lassosum_s-' + para_fix['s'] + '_n.score')
        elif regon == 'va':
            fscore = d / ('lassosum_s-' + para_fix['s'] + '_n.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    else:
        raise NotImplementedError('Unknown method', method)

    return fscore


def sort_fscores_nsnvs(fscores, nsnvs):
    zip_sorted = sorted(zip(fscores, nsnvs), key=lambda pair: int(pair[1]))
    fs = [x for (x, _) in zip_sorted]
    ns = [y for (_, y) in zip_sorted]
    return fs, ns


def load_score_str(f):
    score = pd.read_csv(f, delim_whitespace=True, dtype=str)
    if 'fid' in score:
        score = score.drop(columns=['fid'])
    if 'status' in score:
        score = score.drop(columns=['status'])
    if 'phe' in score:
        score = score.drop(columns=['phe'])
    score = score.rename(columns={'iid': 'id', 'rs': 'score'})
    assert list(score.columns) == ['id', 'score']
    return score


def concat_score(fs, nsnvs, para_concat):
    assert para_concat == 'n'
    scores = []
    for f in fs:
        score = load_score_str(f)
        scores.append(score)

    score_concat = scores[0].loc[:, ['id']].set_index('id').copy()
    # score_concat = scores[0].loc[:, ['id', 'status']].set_index('id').copy()
    for score, n in zip(scores, nsnvs):
        score = score.set_index('id')
        # if not (score_cocnat.set_index('id')['status']==score.set_index('id')['status']).all():
        # if not (score_concat['status'] == score['status']).all():
        # raise RuntimeError('phe is not consistent: n', n)

        # drop 'status'
        score = score[['score']].rename(columns={'score': 'score_n-' + str(n)})

        score_concat = score_concat.join(score)
    logger.debug('score_concat\n{}'.format(score_concat.head()))
    score_concat = score_concat.reset_index(drop=False)
    return score_concat


# this should be faster and more accurate than using regex
def parse_para(f, para):
    # assume suffix only has one '.'
    # ok: '.score_withcov'
    # ng: '.score.withcov'
    base = pathlib.Path(f).stem
    pnames = base.split('_')

    parad = {}
    for pname in pnames:
        if '-' not in pname:
            continue
        ps = pname.split('-')

        if ps[0] in para:
            parad[ps[0]] = ps[1]
    return parad


def uniq_dict(ds):
    uniq = []
    for x in ds:
        if x not in uniq:
            uniq.append(x)
    return uniq


def load_files_withoutn(d, method, regon):
    if method == 'boosting':
        g_fscore = None
    elif method == 'snpnet':
        g_fscore = None
    elif method == 'ldpred':
        if regon is None:
            g_fscore = d / ('ldpred_rho-*.score')
        elif regon == 'va':
            g_fscore = d / ('ldpred_rho-*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'clump':
        if regon is None:
            # exclude '_n-' later
            g_fscore = d / ('clump_p-*_r2-*.score')
        elif regon == 'va':
            g_fscore = d / ('clump_p-*_r2-*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'prscs':
        if regon is None:
            g_fscore = d / ('prscs_a-*_b-*_phi-*.score')
        elif regon == 'va':
            g_fscore = d / ('prscs_a-*_b-*_phi-*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'sbayesr':
        if regon is None:
            g_fscore = d / ('sbayesr.score')
        elif regon == 'va':
            g_fscore = d / ('sbayesr.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'lassosum':
        g_fscore = None
    else:
        raise NotImplementedError('Unknown method', method)

    if g_fscore is None:
        return []
    # do not contain _
    # -> ng
    # g_fscore = d / ('clump_p-*_r2-[^_]+.score')
    # g_fscore = d / ('clump_p-*_r2-[!_]*.score')
    fscores = glob.glob(str(g_fscore))

    fscores = [pathlib.Path(f) for f in fscores]
    # exclude '_n-' here
    fscores = [f for f in fscores if '_n-' not in f.name]

    return fscores


def load_parads_withn(d, method, regon):
    if method == 'boosting':
        paras = []
        g_fscore = d / (method + '_n*.score')
    elif method == 'snpnet':
        paras = []
        if regon is None:
            g_fscore = d / (method + '_n-*.score')
        elif regon == 'tr':
            g_fscore = d / (method + '_n-*.scorecov')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'ldpred':
        g_fscore = None
        # paras = ['rho']
        # g_fscore = d / ('ldpred_rho-*_n*.score')
    elif method == 'clump':
        paras = ['p', 'r2']
        if regon is None:
            g_fscore = d / ('clump_p-*_r2-*_n*.score')
        elif regon == 'va':
            g_fscore = d / ('clump_p-*_r2-*_n*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    elif method == 'prscs':
        g_fscore = None
    elif method == 'sbayesr':
        g_fscore = None
    elif method == 'lassosum':
        paras = ['s']
        if regon is None:
            g_fscore = d / ('lassosum_s-*_n-*.score')
        elif regon == 'va':
            g_fscore = d / ('lassosum_s-*_n-*.scorecov_regonva')
        else:
            raise NotImplementedError('Unknown regon', regon)
    else:
        raise NotImplementedError('Unknown method', method)

    if g_fscore is None:
        return []
    fscores = glob.glob(str(g_fscore))

    logger.debug("fscores {}".format(fscores))

    parads = [parse_para(f, paras) for f in fscores]
    # print('parads', parads)
    # ng for list of dict
    # parads = list(set(parads))
    parads = uniq_dict(parads)

    # print('parads', parads)

    return parads


def concat_paras(dout, dscore, method, para_concat, regon):
    if para_concat != 'n':
        raise NotImplementedError

    logger.debug('Start withoutn')
    # without n; just for changing format of score
    fscores = load_files_withoutn(dscore, method, regon)
    for fscore in fscores:
        logger.debug('fscore: {}'.format(fscore))
        score = load_score_str(fscore)
        fscore_new = dout / fscore.name
        logger.debug('fscore_new:{}'.format(fscore_new))
        # fscore_concat = get_fscore_concat(dout, para_fix={'p': p, 'r2': r2}, para_concat=para_concat)
        score.to_csv(fscore_new, sep='\t', index=False)

    logger.debug('Start withn')
    parads = load_parads_withn(dscore, method, regon)
    # parads = load_parads_withn(dscore, paras=['p', 'r2'])
    # parad={'p':0.01, 'r2':0.1}
    for parad in parads:
        #print('parad', parad)
        fscores = []
        fscores, nsnvs = get_fscores_n(dscore, para_fix=parad, para=para_concat, method=method, regon=regon)
        # fscores, nsnvs = get_fscores(dscore, para_fix={'p': p, 'r2': r2}, para=para_concat)
        #print('nsnvs', nsnvs)
        if len(fscores) == 0:
            continue
        fscores, nsnvs = sort_fscores_nsnvs(fscores, nsnvs)
        score_concat = concat_score(fscores, nsnvs, para_concat=para_concat)
        fscore_concat = get_fscore_concat(dout, para_fix=parad, para_concat=para_concat, method=method, regon=regon)
        logger.debug('fscore_concat: {}'.format(fscore_concat))
        # fscore_concat = get_fscore_concat(dout, para_fix={'p': p, 'r2': r2}, para_concat=para_concat)
        score_concat.to_csv(fscore_concat, sep='\t', index=False)

        # raise NotImplementedError


def main():
    args = argument()

    if args.concat_paras:
        # logger.info('join_snv_hm3]')
        concat_paras(
            pathlib.Path(args.dout),
            pathlib.Path(args.dout_score),
            args.method, args.para_concat,
            args.regon)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
