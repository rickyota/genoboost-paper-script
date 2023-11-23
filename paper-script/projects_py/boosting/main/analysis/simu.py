
import os
import pathlib
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from logging import getLogger

from ...io import filename as iof
from ...wgt import wgt as wgtop
from .plot import genetic_model as plot_genetic_model
from .. import paper

from ....system.logger import logger_setting


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--genetic-model',
                        action='store_true',
                        help='')

    parser.add_argument('--dout',
                        help='output file')

    parser.add_argument('--ddata',
                        help='dir of simu')

    parser.add_argument('--dresult',
                        help='')

    parser.add_argument('--model',
                        help='')

    parser.add_argument('--simu_type',
                        help='')

    # parser.add_argument('--dassoc',
    #                    default=None,
    #                    help='output file')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


# merge to wgtop.
def add_genetic_model(wgt, genetic_threshold_kind='1'):
    def add_relative_score(wgt):
        wgt.loc[:, 'score1_rel'] = wgt['score1'] - wgt['score0']
        wgt.loc[:, 'score2_rel'] = wgt['score2'] - wgt['score0']
        return wgt

    wgt = add_relative_score(wgt)
    # add rec or dom
    wgt = wgtop.estimate_model_score_boost(wgt, score1='score1_rel', score2='score2_rel', genetic_threshold_kind=genetic_threshold_kind)

    # logger.debug('{}'.format(wgt))

    uniq, counts = np.unique(wgt['model_est'].dropna(), return_counts=True)
    d = dict(zip(uniq, counts))

    logger.debug('{}'.format(d))
    return wgt


def format_true_wgt(wgt):
    cols_unchange = 'id,chrom,pos,a1,a2'.split(',')
    wgt = wgt.rename(columns=lambda x: 'true_' + x if x not in cols_unchange else x)
    # cols_drop = 'chrom,pos,a1,a2'.split(',')
    # wgt = wgt.drop(columns=cols_drop)
    # wgt = wgt.rename(columns=lambda x: x + '_true' if x != 'id' else x)
    return wgt


def create_ss_wgt_close_to_true(ss_wgt: pd.DataFrame):

    # extract true-only variants
    # if the variant is not selected by Boosting, use the nearest var in +-1Mbp as the corresponding var.
    true_only = ss_wgt.loc[~ss_wgt['true_score0'].isna(), :].copy()

    logger.debug('true&boosting selected: {}'.format(len(true_only.loc[~true_only['boost_iteration'].isna(), :])))
    logger.debug('true&boosting selected: {}'.format(true_only.loc[~true_only['boost_iteration'].isna(), :]))

    # index should be orrespondent to ss_wgt
    for index in true_only.index:
        if not np.isnan(true_only.loc[index, 'boost_iteration']):
            # both true and boost exists
            true_only.loc[index, 'boost_model_est_close'] = true_only.loc[index, 'boost_model_est']
        else:
            # select closest var
            # print('true-only', true_only.loc[index, :])
            chrom = true_only.loc[index, 'chrom']
            pos = true_only.loc[index, 'pos']
            # print('chom,pos', chrom, pos)
            ss_wgt_near = ss_wgt.loc[ss_wgt['chrom'] == chrom, :]
            # only use < 1Mb
            pos_range = 1000000
            ss_wgt_near = ss_wgt_near.loc[((pos - pos_range) < ss_wgt['pos']) & (ss_wgt['pos'] < (pos + pos_range)), :]

            if len(ss_wgt_near) == 0:
                # no var within +- 1Mb
                # do nth
                continue

            index_i_boost_close = (ss_wgt_near.loc[:, 'pos'] - pos).argmin()
            # https://pandas.pydata.org/docs/user_guide/indexing.html#combining-positional-and-label-based-indexing
            index_boost_close = ss_wgt_near.index[index_i_boost_close]
            # logger.debug('pos : {} vs {}'.format(pos, ss_wgt_near.loc[index_boost_close, 'pos']))

            # assign nearest var
            true_only.loc[index, 'boost_model_est_close'] = ss_wgt_near.loc[index_boost_close, 'boost_model_est']

    return true_only


def create_heatmap_table(ss_wgt, xy):

    data = create_ss_wgt_close_to_true(ss_wgt)

    gmodels = ['overdom', 'dom', 'add', 'rec', 'overrec']

    if xy == ('gwas', 'boost'):
        x_col = 'model_est_gwas'
        y_col = 'boost_model_est_close'
    elif xy == ('true', 'boost'):
        x_col = 'true_model_est'
        y_col = 'boost_model_est_close'
    elif xy == ('true', 'gwas'):
        x_col = 'true_model_est'
        y_col = 'model_est_gwas'
    else:
        raise RuntimeError

    data = data[[y_col, x_col]].copy()

    # data = ss_wgt[[boost_model_est, gwas_model_est]].copy()
    # data = ss_wgt[['boost_model_est', 'model_est_gwas']].copy()
    data.loc[:, 'dummy'] = 1
    # print('data', data.head())

    stat_heatmap = pd.pivot_table(data,
                                  values='dummy',
                                  index=y_col, columns=x_col,
                                  # index='boost_model_est', columns='model_est_gwas',
                                  aggfunc='count')

    logger.debug('data pivot: {}'.format(data))

    def fill_col_index(data, gmodels):
        for gmodel in gmodels:
            # fill index
            if gmodel not in data.index:
                data.loc[gmodel, :] = np.nan
            # fill col
            if gmodel not in data:
                data.loc[:, gmodel] = np.nan
        return data
    stat_heatmap = fill_col_index(stat_heatmap, gmodels)

    stat_heatmap = stat_heatmap.loc[gmodels, gmodels]

    # FIXME: use categorical type for boost_model_est
    # then, these are not necessary
    stat_heatmap = stat_heatmap.fillna(0.0)
    stat_heatmap = stat_heatmap.astype(int)
    logger.debug('stat_heatmap pivot: {}'.format(stat_heatmap))
    return stat_heatmap


def genetic_modeld2():
    return {'overrec': 'Overrec', 'rec': 'Rec', 'add': 'Add', 'dom': 'Dom', 'overdom': 'Overdom'}


def palette_hls_2_genetic_model():

    cb = sns.color_palette("colorblind", 10)
    paletted = {'overdom': cb[6],
                'dom': cb[4],
                'add': cb[2],
                'rec': cb[0],
                'overrec': cb[9]
                }
    return paletted


def validate_heatmap(droot, stat, simu_type, random_name):
    os.makedirs(droot, exist_ok=True)

    with sns.plotting_context('poster', font_scale=1.0):
        xy = ('true', 'boost')
        # true_wgt vs boosting wgt
        heatmap_true = create_heatmap_table(stat, xy)
        fig = plot_genetic_model.plot_genetic_model_gwas_heatmap(heatmap_true,
                                                                 xy,
                                                                 # gmodels=gmodels,
                                                                 gmodeld=genetic_modeld2(),
                                                                 palette=palette_hls_2_genetic_model(),
                                                                 title=simu_type + '_' + random_name.split('_')[1]
                                                                 )

        fname = droot / 'heatmap_true_boost'
        logger.debug('fname {}'.format(fname))
        fig.savefig(fname.parent / (fname.name + '.png'))
        # pio.savefig_png_svg(fig, fname, dpi=1200)
        plt.close()

    with sns.plotting_context('poster', font_scale=1.0):
        xy = ('true', 'gwas')
        heatmap_true = create_heatmap_table(stat, xy)
        fig = plot_genetic_model.plot_genetic_model_gwas_heatmap(heatmap_true,
                                                                 xy,
                                                                 # gmodels=gmodels,
                                                                 gmodeld=genetic_modeld2(),
                                                                 palette=palette_hls_2_genetic_model(),
                                                                 title=simu_type + '_' + random_name.split('_')[1]
                                                                 )

        fname = droot / 'heatmap_true_gwas'
        logger.debug('fname {}'.format(fname))
        fig.savefig(fname.parent / (fname.name + '.png'))
        # pio.savefig_png_svg(fig, fname, dpi=1200)
        plt.close()


def genetic_model_from_wgt(wgt, col):
    # models = gmodels()
    gmodels = ['overdom', 'dom', 'add', 'rec', 'overrec']

    uniq, counts = np.unique(wgt[col].dropna(), return_counts=True)
    # uniq, counts = np.unique(wgt['boost_model_est'].dropna(), return_counts=True)
    d = dict(zip(uniq, counts))

    for k in gmodels:
        if k not in d:
            d[k] = 0

    return d


def stat_genetic_model_from_wgt(wgt_ss, cols):
    # col_true = 'true_model_est'
    # col_gwas = 'model_est_gwas'
    # col_boost = 'boost_model_est_close'

    stats = []
    for col in cols:
        d = genetic_model_from_wgt(wgt_ss, col)
        d |= {'method': col}
        stats.append(d)

    # remove None
    stats = [x for x in stats if x is not None]
    stats = pd.DataFrame.from_records(stats)

    stats = stats.set_index('method')

    return stats


def validate_genetic_model(droot, wgt_ss, simu_type, random_name):
    os.makedirs(droot, exist_ok=True)

    with sns.plotting_context('poster'):
        # with sns.color_palette(mypalette4r, 5):
        cols = ['true_model_est', 'boost_model_est', 'model_est_gwas']
        # true_wgt vs boosting wgt
        stats = stat_genetic_model_from_wgt(wgt_ss, cols)
        print('stats', stats)
        fig = plot_genetic_model.plot_genetic_model(stats,
                                                    # cols,
                                                    gmodeld=genetic_modeld2(),
                                                    title=simu_type + '_' + random_name.split('_')[1]
                                                    )

        fig.tight_layout()
        fname = droot / 'genetic_model_true_boost'
        logger.debug('fname {}'.format(fname))
        fig.savefig(fname.parent / (fname.name + '.png'))
        # pio.savefig_png_svg(fig, fname, dpi=1200)
        plt.close()


def merge_wgt(true_wgt, wgt):
    # cols_drop = 'chrom,pos,a1,a2'.split(',')
    # true_wgt = true_wgt.drop(columns=cols_drop)

    # wgt_ss = pd.merge(true_wgt, wgt, on='id', how='outer')

    cols_merge = 'id,chrom,pos,a1,a2'.split(',')
    wgt_ss = pd.merge(true_wgt, wgt, on=cols_merge, how='outer')
    return wgt_ss


def simu_genetic_model(dout, ddata_simu, dresult, model, simu_type):

    # gmodels = ['add', 'dom', 'rec', 'hetonly']

    cvi = '0'
    kind = 'simu'
    method = 'boosting'

    # simu_types = os.listdir(ddata_simu)

    # TMP
    # simu_types = ['h2add-0.08_h2dom-0.02_ncausal-100']
    # simu_types = ['h2add-0.08_h2dom-0.02_ncausal-1000']
    # simu_types = ['h2add-0.09_h2dom-0.01_ncausal-1000']
    # simu_types = ['h2add-0.45_h2dom-0.05_ncausal-1000']

    # for simu_type in simu_types:
    logger.debug('model: {}'.format(simu_type))
    dwgt = ddata_simu / simu_type / 'true-wgt'
    for fwgt_base in os.listdir(dwgt):
        random_name = fwgt_base.split('.wgt')[0]

        #if random_name != 'random_0':
        #    print("TMP")
        #    continue

        fwgt = dwgt / fwgt_base
        true_wgt = pd.read_csv(fwgt, sep='\t')
        true_wgt = add_genetic_model(true_wgt, '1')
        true_wgt = format_true_wgt(true_wgt)

        # tmp
        uniq, counts = np.unique(true_wgt['true_model_est'].dropna(), return_counts=True)
        d = dict(zip(uniq, counts))
        logger.debug('true wgt {}'.format(d))

        fwgt = iof.file_wgt_boosting_best_para(dresult, method, random_name, cvi, [kind, simu_type, model])
        #  best para
        wgt = wgtop.load_wgt_method(fwgt, method)
        wgt = wgtop.format_wgt_boost(wgt)
        print('wgt', wgt)

        # tmp
        uniq, counts = np.unique(wgt['boost_model_est'].dropna(), return_counts=True)
        d = dict(zip(uniq, counts))
        logger.debug('wgt {}'.format(d))

        wgt_ss = merge_wgt(true_wgt, wgt)

        # use 'update_ssd_wgt' instead
        ssd = paper.load_ssd_raw(dresult, random_name, cvi, mid_path=[kind, simu_type])
        wgt_ss = paper.merge_wgt_ssd(wgt_ss, ssd, 'left')
        wgt_ss = paper.judge_genetic_model_gwas(wgt_ss)

        print('wgt_ss', wgt_ss)
        print('wgt_ss', wgt_ss.columns)

        droot = dout / simu_type / random_name
        validate_heatmap(droot, wgt_ss, simu_type, random_name)

        validate_genetic_model(droot, wgt_ss, simu_type, random_name)

        # TMP
        raise RuntimeError

    # logger.debug('Done!')


def main():
    args = argument()

    if args.genetic_model:
        simu_genetic_model(pathlib.Path(args.dout),
                           pathlib.Path(args.ddata),
                           pathlib.Path(args.dresult),
                           args.model, args.simu_type)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
