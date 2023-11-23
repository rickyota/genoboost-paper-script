
import argparse
from logging import getLogger
import pandas as pd
import numpy as np


import statsmodels.api as sm

import seaborn as sns
from matplotlib import pyplot as plt

import matplotlib
import matplotlib.ticker as ticker
import seaborn as sns
matplotlib.use('Agg')
# plt.rcParams['patch.edgecolor'] = 'none'
# plt.rcParams['lines.markeredgewidth'] = 0
# this might reset configs
# sns.set(style="ticks")
sns.set_style("ticks")

from ....system.logger import logger_setting


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--dout',
                        default=None,
                        help='output file')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def linreg_quan_vs_quan(y, X):
    # ys = (1, -1)
    # ycase, ycont = ys
    # y = [ycase] * table[0, 0] + [ycont] * table[0, 1] + [ycase] * table[1, 0] + [ycont] * table[1, 1]
    # X = [1] * table[0, 0] + [1] * table[0, 1] + [0] * table[1, 0] + [0] * table[1, 1]

    y = np.array(y)
    X = pd.DataFrame(X, columns=['X'])
    X.loc[:, 'const'] = 1

    # print(y)
    # print(X)
    # X = sm.add_constant(X)
    # print("add const", X.shape)

    logit = sm.OLS(y, X).fit()
    # logit = sm.Logit(y, X).fit()
    logger.debug('logit.summary: {}'.format(logit.summary()))

    # print("coeff",logit.params.iloc[0])
    alpha = logit.params.iloc[0]
    const = logit.params.iloc[1]  # const is the last

    se_alpha = logit.bse.iloc[0]
    se_const = logit.bse.iloc[1]

    pval_alpha = logit.pvalues.iloc[0]

    return alpha, const, se_alpha, se_const, pval_alpha


def count(y, X, y_label, x):
    # print('y', (y == y_label))
    filt = (np.array(y) == y_label) & (np.array(X) == x)
    # filt = np.logical_and((y == y_label), (X == x))
    # print('filt', filt)
    return filt.sum()


def prepare_samples():
    samples_n = 100
    # should use y=(1,0) since to comapre with ss_logistic
    # binomial
    # x=(0,30)
    mean = 1.0
    std = 1.0
    x = np.random.normal(mean, std, size=samples_n)
    # print(x)
    y = - x
    # max=1
    # (x, y) = (x / n - 0.5, y / n - 0.5)
    # case, control
    label = [1] * samples_n + [0] * samples_n
    X = list(x) + list(y)

    # 1 if case and top10% cov, or control and  bottom 10%?
    case_prop = 0.1
    x_topn = int((samples_n) * (1 - case_prop))
    x_top = np.sort(x)[x_topn]
    logger.debug('x max min: {}, {}'.format(np.max(x), np.min(x)))
    logger.debug('x_top: {}'.format(x_top))
    X_snv_case = [1 if x_ >= x_top else 0 for x_ in x]

    # pattern1: both case and cont have snv
    # y_top = np.sort(y)[-1]
    # X_snv_cont = [1 if y_ >= y_top else 0 for y_ in y]
    # X_snv = list(X_snv_case) + list(X_snv_cont)

    # pattern 1.5: more than one cont have snv
    # y_topn = int((samples_n) * case_prop)
    # y_top = np.sort(y)[y_topn]

    # pattern2: only case have snv
    X_snv = list(X_snv_case) + [0] * samples_n

    # print('X', X)
    # print('X_snv', X_snv)

    count_table_1 = [count(label, X_snv, 1, 1), count(label, X_snv, 0, 1)]
    count_table_0 = [count(label, X_snv, 1, 0), count(label, X_snv, 0, 0)]
    logger.debug('count snv case/cont: {}, {}'.format(count_table_1[0], count_table_1[1]))
    logger.debug('count no snv case/cont: {}, {}'.format(count_table_0[0], count_table_0[1]))

    return label, X, X_snv


"""
# ok param to save
def prepare_samples():
    samples_n = 1000
    # should use y=(1,0) since to comapre with ss_logistic
    # binomial
    # x=(0,30)
    n = 30
    prop = 0.55
    x = np.random.binomial(n, prop, samples_n)
    # print(x)
    y = n - x
    # max=1
    (x, y) = (x / n - 0.5, y / n - 0.5)
    # case, control
    label = [1] * samples_n + [0] * samples_n
    X = list(x) + list(y)

    # 1 if case and top10% cov, or control and  bottom 10%?
    case_prop = 0.05
    x_topn = int((samples_n) * (1 - case_prop))
    x_top = np.sort(x)[x_topn]
    print('x max min', np.max(x), np.min(x))
    print('x_top', x_top)
    X_snv_case = [1 if x_ >= x_top else 0 for x_ in x]

    # pattern1: both case and cont have snv
    #y_top = np.sort(y)[-1]
    #X_snv_cont = [1 if y_ >= y_top else 0 for y_ in y]
    #X_snv = list(X_snv_case) + list(X_snv_cont)

    # pattern 1.5: more than one cont have snv
    ##y_topn = int((samples_n) * case_prop)
    ##y_top = np.sort(y)[y_topn]

    # pattern2: only case have snv
    X_snv = list(X_snv_case) + [0] * samples_n

    #print('X', X)
    #print('X_snv', X_snv)

    count_table_1 = [count(label, X_snv, 1, 1), count(label, X_snv, 0, 1)]
    count_table_0 = [count(label, X_snv, 1, 0), count(label, X_snv, 0, 0)]
    print('count snv case/cont', count_table_1[0], count_table_1[1])
    print('count no snv case/cont', count_table_0[0], count_table_0[1])

    return label, X, X_snv
"""


def logreg(y, X):
    # N = np.sum(table)
    # y = np.empty_like(N)
    # X = np.empty_like(N)
    # y = [1] * table[0, 0] + [0] * table[0, 1] + [1] * table[1, 0] + [0] * table[1, 1]
    # X = [1] * table[0, 0] + [1] * table[0, 1] + [0] * table[1, 0] + [0] * table[1, 1]
    y = np.array(y)
    X = pd.DataFrame(X, columns=['X'])
    X.loc[:, 'const'] = 1
    # print(y)
    # print(X)
    # X = sm.add_constant(X)
    # print("add const", X.shape)

    logit = sm.Logit(y, X).fit()
    logger.debug('logit.summary: {}'.format(logit.summary()))

    # print("coeff",logit.params.iloc[0])
    alpha = logit.params.iloc[0]
    const = logit.params.iloc[1]  # const is the last

    se_alpha = logit.bse.iloc[0]
    se_const = logit.bse.iloc[1]

    pval_alpha = logit.pvalues.iloc[0]

    return alpha, const, se_alpha, se_const, pval_alpha


def logreg_cov(y, X, cov):
    # X, y, cov = create_X_y_cov(table, cov_case, cov_cont)
    # y=[1]*table[0,0]+[0]*table[0,1]+[1]*table[1,0]+[0]*table[1,1]+[1]*table[2,0]+[0]*table[2,1]
    # X=[2]*table[0,0]+[2]*table[0,1]+[1]*table[1,0]+[1]*table[1,1]+[0]*table[2,0]+[0]*table[2,1]

    y = np.array(y)
    X = pd.DataFrame(np.array([X, cov]).T, columns=['X', 'cov'])
    X.loc[:, 'const'] = 1

    # print(y)
    # print(X)
    # X = sm.add_constant(X)
    # print("add const", X.shape)

    logit = sm.Logit(y, X).fit()
    logger.debug('logit.summary: {}'.format(logit.summary()))
    #print(logit.summary())

    # print("coeff",logit.params.iloc[0])
    alpha = logit.params.iloc[0]
    alpha_cov = logit.params.iloc[1]
    const = logit.params.iloc[2]  # const is the last

    se_alpha = logit.bse.iloc[0]
    se_alpha_cov = logit.bse.iloc[1]
    se_const = logit.bse.iloc[2]

    pval_alpha = logit.pvalues.iloc[0]

    return alpha, alpha_cov, const, se_alpha, se_alpha_cov, se_const, pval_alpha


def w_sum(ws, y, X, y_label, x):
    # print('y', (y == y_label))
    filt = (y == y_label) & (X == x)
    # filt = np.logical_and((y == y_label), (X == x))
    # print('filt', filt)
    return ws[filt].sum()


def run_modelfree(y, X, cov):
    y = np.array(y)
    X = np.array(X)
    cov = np.array(cov)
    print('logreg for cov')
    alpha_cov, const_cov, _, _, _ = logreg(y, cov)
    print('alpha_cov, const_cov', alpha_cov, const_cov)

    score = np.array([alpha_cov * x + const_cov for x in cov])
    ws = 1 / (1 + np.exp(y * score))
    ws = ws / ws.sum()
    print('ws', ws[0:100])
    # print('y', y)
    # print('X', X)

    w_table_1 = [w_sum(ws, y, X, y_label=1, x=1), w_sum(ws, y, X, 0, 1)]
    w_table_0 = [w_sum(ws, y, X, 1, 0), w_sum(ws, y, X, 0, 0)]
    if w_table_1[1] == 0:
        print('eps')
        # case/2
        ws_med = np.median(ws[y == 1]) / 2
        print('eps=', ws_med)
        w_table_1[1] = ws_med

    print('w_table_1', w_table_1)
    print('w_table_0', w_table_0)
    score_1 = np.log(w_table_1[0] / w_table_1[1]) / 2
    score_0 = np.log(w_table_0[0] / w_table_0[1]) / 2
    alpha = score_1 - score_0
    print('modelfree alpha', alpha)

    if alpha >= 0:
        print("Great modelfree positive")
    else:
        print("really...? modelfree negative")

    return alpha


def plot(fout_plot, X, label, X_snv, alpha, const, alpha_snv, const_snv, res):
    figsize = (8, 8)

    s = 60
    # s = 10
    alpha_c = 0.5
    lw = 2
    m0 = 'o'  # case, no snv
    m1 = 'x'  # 'x'  # case, with snv
    m2 = '+'  # '+'
    c0 = 'red'
    c1 = 'orange'
    c2 = 'blue'

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.scatter(X, label, s=s)
    ax.set_ylim((-0.5, 1.5))
    fig.tight_layout()
    fig.savefig(fout_plot + '.input.png')
    fig.savefig(fout_plot + '.input.pdf')

    # change color for X_snv
    X_have_snv = [x for (x_snv, x) in zip(X_snv, X) if x_snv == 1]
    label_have_snv = [l_ for (x_snv, l_) in zip(X_snv, label) if x_snv == 1]
    X_case = [x for (x, x_snv, lab) in zip(X, X_snv, label) if (lab == 1) and (x_snv == 0)]
    label_case = [lab for (x_snv, lab) in zip(X_snv, label) if (lab == 1) and (x_snv == 0)]
    # X_case = [x for (x, lab) in zip(X, label) if lab == 1]
    # label_case = [lab for lab in label if lab == 1]
    X_cont = [x for (x, lab) in zip(X, label) if lab == 0]
    label_cont = [lab for lab in label if lab == 0]
    print('label_have_snv', label_have_snv)

    def plot_dot(ax, x, y, dtype=0):
        if dtype == 0:
            c, m = c0, m0
            label = 'Case without an SNV'
            ax.scatter(x, y, s=s, marker=m, alpha=alpha_c, lw=lw, facecolors='None', edgecolors=c, label=label)
        elif dtype == 1:
            c, m = c1, m1
            label = 'Case with an SNV'
            ax.scatter(x, y, s=s, c=c, marker=m, label=label)
        elif dtype == 2:
            c, m = c2, m2
            label = 'Control without an SNV'
            ax.scatter(x, y, s=s, c=c, marker=m, alpha=alpha_c, label=label)
        # ax.scatter(x, y, s=s, marker=m, alpha=alpha, lw=lw, facecolors='None', edgecolors=c)
        # ax.scatter(x, y, s=s, c=c, marker=m, alpha=alpha, lw=lw, edgecolors='white')
        return ax

    # for legend
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax = plot_dot(ax, X_case, label_case, 0)
    ax = plot_dot(ax, X_have_snv, label_have_snv, 1)
    ax = plot_dot(ax, X_cont, label_cont, 2)
    ax.set_ylim((-0.5, 1.5))
    ax.set(ylabel='Label ($y_i$)', xlabel='Covariate ($cov_i$)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(fout_plot + '.input.color.legend.png')
    fig.savefig(fout_plot + '.input.color.legend.pdf')

    # regression
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax = plot_dot(ax, X_case, label_case, 0)
    ax = plot_dot(ax, X_cont, label_cont, 2)
    ax = plot_dot(ax, X_have_snv, label_have_snv, 1)
    # ax.scatter(X_case, label_case, s=s, c=c0, marker=m0, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # ax.scatter(X_cont, label_cont, s=s, c=c2, marker=m2, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # ax.scatter(X_have_snv, label_have_snv, s=s, c=c1, marker=m1, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # for below
    xmin, xmax = ax.get_xlim()
    xs_lin = (xmin, xmax)
    ys_lin = (const + xmin * alpha, const + xmax * alpha)
    ax.plot(xs_lin, ys_lin, c=c0)
    ax.set_xlim((-4.0, 4.0))
    ylims = ax.get_ylim()
    ax.set(ylabel='Label ($y_i$)', xlabel='Covariate ($cov_i$)')
    fig.tight_layout()
    fig.savefig(fout_plot + '.reg1.png')
    fig.savefig(fout_plot + '.reg1.pdf')

    # input plot
    # here since want to make ylim the same as above
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax = plot_dot(ax, X_case, label_case, 0)
    ax = plot_dot(ax, X_cont, label_cont, 2)
    ax = plot_dot(ax, X_have_snv, label_have_snv, 1)
    # ax.scatter(X_case, label_case, s=s, c=c0, marker=m0, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # ax.scatter(X_cont, label_cont, s=s, c=c2, marker=m2, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # ax.scatter(X_have_snv, label_have_snv, s=s, c=c1, marker=m1, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    ax.set_xlim((-4.0, 4.0))
    ax.set_ylim(ylims)
    # ax.set_ylim((-0.5, 1.5))
    ax.set(ylabel='Label ($y_i$)', xlabel='Covariate ($cov_i$)')
    fig.tight_layout()
    fig.savefig(fout_plot + '.input.color.png')
    fig.savefig(fout_plot + '.input.color.pdf')

    # residual
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    res_have = [res_ for (x_snv, res_) in zip(X_snv, res) if x_snv == 1]
    X_have = [1] * len(res_have)

    X_snv_case = [x for (x, lab) in zip(X_snv, label) if (lab == 1) and (x == 0)]
    res_case = [r for (r, x, lab) in zip(res, X_snv, label) if (lab == 1) and (x == 0)]
    # X_snv_case = [x for (x, lab) in zip(X_snv, label) if lab == 1]
    # res_case = [r for (r, lab) in zip(res, label) if lab == 1]
    X_snv_cont = [x for (x, lab) in zip(X_snv, label) if lab == 0]
    res_cont = [r for (r, lab) in zip(res, label) if lab == 0]

    ax = plot_dot(ax, X_snv_case, res_case, 0)
    ax = plot_dot(ax, X_snv_cont, res_cont, 2)
    ax = plot_dot(ax, X_have, res_have, 1)
    # ax.scatter(X_snv_case, res_case, s=s, c=c0, marker=m0, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # ax.scatter(X_snv_cont, res_cont, s=s, c=c2, marker=m2, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    # ax.scatter(X_have, res_have, s=s, c=c1, marker=m1, alpha=alpha, lw=lw, facecolors='None', edgecolors='white')
    xmin, xmax = ax.get_xlim()
    xs_lin = (xmin, xmax)
    ys_lin = (const_snv + xmin * alpha_snv, const_snv + xmax * alpha_snv)
    ax.plot(xs_lin, ys_lin, c=c0)
    ax.set(ylabel='Residual ($res_i$)', xlabel='SNV ($x_i$)')
    fig.tight_layout()
    fig.savefig(fout_plot + '.reg2.png')
    fig.savefig(fout_plot + '.reg2.pdf')


def create_example_l2boost_fail(dout):

    fout_plot = dout + "plot"

    np.random.seed(53)

    label, X, X_snv = prepare_samples()
    # samples_n = 1000
    # label, X, X_snv = prepare_samples(samples_n)

    alpha, const, _, _, _ = linreg_quan_vs_quan(label, X)
    print(alpha, const)

    score = [alpha * x + const for x in X]
    res = [label_ - score_ for (label_, score_) in zip(label, score)]

    alpha_snv, const_snv, _, _, _ = linreg_quan_vs_quan(res, X_snv)
    print(alpha_snv, const_snv)

    if alpha_snv > 0:
        print('Oh,.. not ideal since alpha>=0')
    else:
        print('Great since alpha<=0')

    print('\n\nplot')
    plot(fout_plot, X, label, X_snv, alpha, const, alpha_snv, const_snv, res)
    print('\n\n')

    # print('X_snv only')
    # print('X_snv', X_snv)
    # print('label', label)
    # print('X_snv', np.unique(X_snv, return_counts=True))
    # print('label', np.unique(label, return_counts=True))
    # alpha, const, _, _, _ = logreg(label, X_snv)
    # print(alpha, const)

    print('X_snv log reg with cov')
    alpha, alpha_cov, const, _, _, _, _ = logreg_cov(label, X_snv, X)
    print(alpha, alpha_cov, const)

    if alpha < 0:
        print('Why negative?')
    else:
        print('Great positive')

    run_modelfree(label, X_snv, X)

    print('Done!')


def main():
    args = argument()

    create_example_l2boost_fail(args.dout)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
