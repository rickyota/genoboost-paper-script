
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
from logging import getLogger

# from ..validate.paper import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)


def _plot_manhattan_colorbar(ss, y,
                             significant=None,
                             chrom=None,
                             ax=None,
                             gradcolor=None,
                             xaxis_skip=[],
                             ylabel=None,
                             xlim=None,
                             ylim=None,  # somehow, ylim=(0,None) is wierd
                             # yneg=False,
                             # ymax=None,
                             ):
    """
    if gradcolor=True, use 'boost_iteration_sparse'. if the column does not exist, create one.

    if you want to indicate colormap on ex. iteration, indicate 'gradcolor=True', for others indicated 'gradcolor=serial'
    """

    # logger.debug('ylim {}'.format(ylim))

    im = None

    if chrom is None:
        if ss[y].isnull().all():
            return ax, im
    else:
        if ss.loc[ss['chrom'] == chrom, y].isnull().all():
            return ax, im

    if xlim is True:
        raise RuntimeError('use "exact" instead.')
    if xlim == 'exact':
        # xlim: (0, last), default add empty before and last
        xmin = ss['serial'].min()
        xmax = ss['serial'].max()
        xlim = (xmin, xmax)

    color = ['#000088']
    color_g = ['#b8b8b8']

    x_labels = []
    x_labels_pos = []
    if chrom is None:
        for chrom in range(1, 23):
            ss_chrom = ss[ss['chrom'] == chrom]
            if chrom in xaxis_skip:
                x_labels.append('')
            else:
                x_labels.append(chrom)
            x_labels_pos.append(ss_chrom['serial'].iloc[0])
            # x_labels_pos.append(ss_chrom['serial'].iloc[0])
        ss_plot = ss.copy()
    else:
        ss_chrom = ss[ss['chrom'] == chrom]
        x_labels.append(chrom)
        x_labels_pos.append(ss_chrom['serial'].iloc[0])
        ss_plot = ss_chrom.copy()

    if significant is None:
        linewidth = 0
        s = 3
        ax = ss_plot.plot(kind='scatter', x='serial', y=y,
                          color=color, ax=ax, s=s,
                          linewidth=linewidth,
                          rasterized=True,
                          marker='.'
                          )
    else:
        if len(ss_plot[ss_plot[y].abs() > significant]) > 0:
            # s = 10
            linewidth = 0
            s = 3

            # s = min(20, 1000000 / len(ss_plot))
            # TODO: now only plot ss[y]>significant
            if gradcolor is None:
                ax = ss_plot[ss_plot[y].abs() > significant].plot(kind='scatter', x='serial', y=y,
                                                                  color=color, ax=ax, s=s,
                                                                  linewidth=linewidth,
                                                                  rasterized=True,
                                                                  marker='.'
                                                                  # marker='o'
                                                                  )
            else:
                raise NotImplementedError
                if gradcolor is True:
                    gradcolor_ = 'boost_iteration_sparse'
                    if gradcolor_ not in ss_plot:
                        ss_plot.loc[:, gradcolor_] = (
                            ss_plot['boost_iteration'] / 100).astype(np.int) * 100
                else:
                    # ex. color = 'boost_iteration'
                    gradcolor_ = gradcolor
                """
                ax = ss_plot[ss_plot[y].abs() > significant].plot(kind='scatter', x='serial', y=y,
                                                                c=gradcolor_, colormap='jet',
                                                                ax=ax, s=100000/len(ss_plot))
                """
                df_sig = ss_plot[ss_plot[y].abs() > significant]
                im = ax.scatter(df_sig['serial'], df_sig[y],
                                c=df_sig[gradcolor_], cmap='jet', s=s,
                                linewidth=linewidth,
                                rasterized=True
                                )

        if len(ss_plot[ss_plot[y].abs() <= significant]) > 0:
            # print("plot under sig")
            s = min(10, 10000000 / len(ss_plot))
            linewidth = 0
            ax = ss_plot[ss_plot[y].abs() <= significant].plot(kind='scatter', x='serial', y=y,
                                                               color=color_g, ax=ax, s=s,
                                                                    linewidth=linewidth,
                                                               rasterized=True
                                                               )

    # print("x_labels_pos", x_labels_pos)
    # print("x_labels", x_labels)
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    if chrom is None:
        # temporary
        # ax.set_xlim([0, len(ss)])
        ax.set_xlim((0, ss['serial'].max()))
    else:
        # ??
        # ax.set_xlim([0, len(group)])
        pass

    # if ymax is None:
    #    if chrom is None:
    #        ymax = ss[y].max(skipna=True) * 1.05
    #        # ymax = ss[y].max(numeric_only=True)*1.05
    #    else:
    #        ymax = ss.loc[ss['chrom'] == chrom, y].max(skipna=True) * 1.05
    #        # ymax = ss.loc[ss['chrom'] == chrom, y].max()*1.05
    # else:
    #    ymax = ymax
    #    print("ymax vs max", ymax, ss[y].max(skipna=True))

    if ylim is not None:
        # logger.debug('ylim {}'.format(ylim))
        if ylim == 'ymax':
            ymax = ss[y].max(skipna=True) * 1.05
            ymin = ss[y].min(skipna=True)
            if ymin > 0.0:
                ymin = ymin * 1.05
            else:
                ymin = ymin * 0.95
            ax.set_ylim([ymin, ymax])
        elif ylim == 'ymax_pos':
            ymax = ss[y].max(skipna=True) * 1.05
            ax.set_ylim([0, ymax])
        else:
            logger.debug('ylim {}'.format(ylim))
            raise NotImplementedError('ylim', ylim)

        # error: ylim=(0,None)
        # ax.set_ylim(ylim)

    # if not yneg:
    #    if ylim is not None:
    #        ax.set_ylim([0, ylim])
    #    else:
    #        ax.set_ylim([0, ymax])
    # else:
    #    if ylim is not None:
    #        ax.set_ylim([-ylim, ylim])
    #    else:
    #        ax.set_ylim([-ymax, ymax])

    if xlim is not None:
        ax.set_xlim(xlim)

    ax.set_xlabel('Chromosome')
    # ax.set_xlabel('Chromosome', fontsize=24)

    # remove minor xticks
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    if ylabel is None:
        ax.set_ylabel(y)
        # ax.set_ylabel(y, fontsize=24)
    else:
        ax.set_ylabel(ylabel)
        # ax.set_ylabel(ylabel, fontsize=24)

    if significant is not None:
        ax.axhline(y=significant, color='gray', linestyle='-', linewidth=2.0)

    # ax.set_xticks(fontsize=32, rotation=60)
    # ax.set_yticks(fontsize=32)
    # plt.axhline(y=significant, color='gray', linestyle='-', linewidth = 2.0)
    # plt.xticks(fontsize=32, rotation=60)
    # plt.yticks(fontsize=32)

    # return fig,ax
    return ax, im


# def plot_manhattan_colorbar_old(ss, y, significant, ymax=None, chrom=None, ax=None, yneg=False, gradcolor=None,
#                            xaxis_skip=[], ylabel=None, xlim=None, ylim=None):
#    """
#    if gradcolor=True, use 'boost_iteration_sparse'. if the column does not exist, create one.
#
#    if you want to indicate colormap on ex. iteration, indicate 'gradcolor=True', for others indicated 'gradcolor=serial'
#    """
#
#    im = None
#
#    if chrom is None:
#        if ss[y].isnull().all():
#            return ax, im
#    else:
#        if ss.loc[ss['chrom'] == chrom, y].isnull().all():
#            return ax, im
#
#    if ymax is None:
#        if chrom is None:
#            ymax = ss[y].max(skipna=True) * 1.05
#            # ymax = ss[y].max(numeric_only=True)*1.05
#        else:
#            ymax = ss.loc[ss['chrom'] == chrom, y].max(skipna=True) * 1.05
#            # ymax = ss.loc[ss['chrom'] == chrom, y].max()*1.05
#    else:
#        ymax = ymax
#        print("ymax vs max", ymax, ss[y].max(skipna=True))
#
#    if xlim is True:
#        # xlim: (0, last), default add empty before and last
#        xmin = ss['serial'].min()
#        xmax = ss['serial'].max()
#        xlim = (xmin, xmax)
#
#    color = ['#000088']
#    color_g = ['#b8b8b8']
#
#    significant = significant
#
#    x_labels = []
#    x_labels_pos = []
#    if chrom is None:
#        for chrom in range(1, 23):
#            ss_chrom = ss[ss['chrom'] == chrom]
#            if chrom in xaxis_skip:
#                x_labels.append('')
#            else:
#                x_labels.append(chrom)
#            x_labels_pos.append(ss_chrom['serial'].iloc[0])
#            # x_labels_pos.append(ss_chrom['serial'].iloc[0])
#        ss_plot = ss.copy()
#    else:
#        ss_chrom = ss[ss['chrom'] == chrom]
#        x_labels.append(chrom)
#        x_labels_pos.append(ss_chrom['serial'].iloc[0])
#        ss_plot = ss_chrom.copy()
#
#    if len(ss_plot[ss_plot[y].abs() > significant]) > 0:
#        # s = 10
#        linewidth = 0
#        s = 3
#
#        # s = min(20, 1000000 / len(ss_plot))
#        # TODO: now only plot ss[y]>significant
#        if gradcolor is None:
#            ax = ss_plot[ss_plot[y].abs() > significant].plot(kind='scatter', x='serial', y=y,
#                                                              color=color, ax=ax, s=s,
#                                                              linewidth=linewidth,
#                                                              rasterized=True,
#                                                              marker='.'
#                                                              # marker='o'
#                                                              )
#        else:
#            if gradcolor is True:
#                gradcolor_ = 'boost_iteration_sparse'
#                if gradcolor_ not in ss_plot:
#                    ss_plot.loc[:, gradcolor_] = (
#                        ss_plot['boost_iteration'] / 100).astype(np.int) * 100
#            else:
#                # ex. color = 'boost_iteration'
#                gradcolor_ = gradcolor
#            """
#            ax = ss_plot[ss_plot[y].abs() > significant].plot(kind='scatter', x='serial', y=y,
#                                                              c=gradcolor_, colormap='jet',
#                                                              ax=ax, s=100000/len(ss_plot))
#            """
#            df_sig = ss_plot[ss_plot[y].abs() > significant]
#            im = ax.scatter(df_sig['serial'], df_sig[y],
#                            c=df_sig[gradcolor_], cmap='jet', s=s,
#                            linewidth=linewidth,
#                            rasterized=True
#                            )
#
#    if len(ss_plot[ss_plot[y].abs() <= significant]) > 0:
#        # print("plot under sig")
#        s = min(10, 10000000 / len(ss_plot))
#        linewidth = 0
#        ax = ss_plot[ss_plot[y].abs() <= significant].plot(kind='scatter', x='serial', y=y,
#                                                           color=color_g, ax=ax, s=s,
#                                                                linewidth=linewidth,
#                                                           rasterized=True
#                                                           )
#
#    # print("x_labels_pos", x_labels_pos)
#    # print("x_labels", x_labels)
#    ax.set_xticks(x_labels_pos)
#    ax.set_xticklabels(x_labels)
#    if chrom is None:
#        # temporary
#        # ax.set_xlim([0, len(ss)])
#        ax.set_xlim((0, ss['serial'].max()))
#    else:
#        # ??
#        # ax.set_xlim([0, len(group)])
#        pass
#
#    if not yneg:
#        if ylim is not None:
#            ax.set_ylim([0, ylim])
#        else:
#            ax.set_ylim([0, ymax])
#    else:
#        if ylim is not None:
#            ax.set_ylim([-ylim, ylim])
#        else:
#            ax.set_ylim([-ymax, ymax])
#
#    if xlim is not None:
#        ax.set_xlim(xlim)
#
#    ax.set_xlabel('Chromosome')
#    # ax.set_xlabel('Chromosome', fontsize=24)
#
#    # remove minor xticks
#    ax.xaxis.set_minor_locator(ticker.NullLocator())
#
#    if ylabel is None:
#        ax.set_ylabel(y)
#        # ax.set_ylabel(y, fontsize=24)
#    else:
#        ax.set_ylabel(ylabel)
#        # ax.set_ylabel(ylabel, fontsize=24)
#    ax.axhline(y=significant, color='gray', linestyle='-', linewidth=2.0)
#
#    # ax.set_xticks(fontsize=32, rotation=60)
#    # ax.set_yticks(fontsize=32)
#    # plt.axhline(y=significant, color='gray', linestyle='-', linewidth = 2.0)
#    # plt.xticks(fontsize=32, rotation=60)
#    # plt.yticks(fontsize=32)
#
#    # return fig,ax
#    return ax, im


def plot_manhattan(
        ss,
        figtype=None,
        methods=None,  # to limit methods and order
        methodd=None,
        title=None,
        # no_x=False,
        ycol_method='effect',  # or 'effect_abs'
        ylim=None,
        markersize=None,
        ylabel=None,
        xaxis_skip=[],  # None?
):

    ssobs_wgt = ss

    # TODO: plot iter as color
    # ssobs_wgt.loc[:, 'boost_iteration_sparse'] = ( ssobs_wgt['boost_iteration'] / 100).astype(np.int) * 100

    ssobs_wgt.loc[:, 'serial'] = range(len(ssobs_wgt))

    nrows = 1 + len(methods)

    _, fig_aspect = pcommon.parse_figtype(figtype)
    figsize = (fig_aspect['height'], fig_aspect['height'] * fig_aspect['aspect'])

    fig, axs = plt.subplots(nrows=nrows, figsize=figsize, sharex=True)

    # pval
    _, _ = _plot_manhattan_colorbar(ssobs_wgt, y='-log10(P)',
                                    # significant=0.0,
                                    ax=axs[0], xaxis_skip=xaxis_skip,
                                    xlim='exact', ylim=ylim)

    for ax, method in zip(axs[1:], methods):
        y = method + '_' + ycol_method

        if method.startswith('boosting_'):
            ylabel = methodd[method]
            ax, _ = _plot_manhattan_colorbar(ssobs_wgt,
                                             y=y,
                                             # y=method + '_' + ycol_method,
                                             # y=method + '_effect_abs',
                                             # significant=0.0,
                                             ax=ax, ylabel=ylabel, xaxis_skip=xaxis_skip,
                                             xlim='exact', ylim=ylim)
            # significant=0.0, ax=ax, ylabel=methodd['boosting'], xaxis_skip=xaxis_skip, xlim=True)

        else:
            ax, _ = _plot_manhattan_colorbar(ssobs_wgt,
                                             y=y,
                                             # y=method + '_' + ycol_method,
                                             # y=method + '_effect_abs',
                                             # significant=0.0,
                                             ax=ax, ylabel=methodd[method], xaxis_skip=xaxis_skip,
                                             xlim='exact', ylim=ylim)

        # add modelsize
        model_size = len(ssobs_wgt[~ssobs_wgt[y].isna()])
        # model_size = len(ssobs_wgt[~ssobs_wgt[method + '_' + ycol_method].isna()])
        ylim_max = ax.get_ylim()[1]
        # logger.debug('ylim method {}, {}'.format(method, ylim_max))
        ax.text(0.01, ylim_max * 0.95, '#selected variants = ' + str(model_size),
                fontsize='x-small',
                horizontalalignment='left', verticalalignment='top',
                )

    # fig.colorbar(im, ax=[ax1, ax2, ax3, ax4, ax5])
    fig.suptitle(title)

    return fig


def get_figsize_manhattan_plt_overlay_pos_multi(figtype, plot_gene, gm_n, ncols=1, sharex='col', sharey='row'):
    # +1 for gene
    if plot_gene:
        # gm_n = len(genetic_models)
        if figtype == '1':
            figwidth = 8
            figsize = (figwidth * ncols, 5 * gm_n + 2.5)
            # [1.0] for blank padding axis
            height_ratios = [5] * gm_n + [1.0] + [2.5]
        elif figtype == '2':
            figwidth = 5
            figsize = (figwidth * ncols, 2 * gm_n + 1.0)
            height_ratios = [2] * gm_n + [1.0]
        elif figtype == '3':
            raise NotImplementedError('need pad')
            figwidth = 8
            # 1.0 for title
            figsize = (figwidth * ncols, 1.0 + 2.5 * gm_n + 2.5)
            height_ratios = [2.5] * gm_n + [2.5]
        elif figtype == '4':
            raise NotImplementedError('need pad')
            figwidth = 6
            # 1.0 for title
            figsize = (figwidth * ncols, 1.0 + 4 * gm_n + 2.5)
            height_ratios = [2.5] * gm_n + [2.5]
        elif figtype == '5':
            figwidth = 7
            figsize = (figwidth * ncols, 4 * gm_n + 2.5)
            # [1.0] for blank padding axis
            height_ratios = [4] * gm_n + [1.0] + [2.5]
        # elif figtype == '6':
        #    figwidth = 7
        #    figsize = (figwidth * ncols, 4.5 * gm_n + 2.5)
        #    # [1.0] for blank padding axis
        #    height_ratios = [4.5] * gm_n + [1.0] + [2.5]
        else:
            raise RuntimeError()
        # +2 for gene and black padding axis bet. man and gene
        fig, axs = plt.subplots(nrows=gm_n + 2, ncols=ncols, figsize=figsize, sharex=sharex,
                                sharey=sharey, gridspec_kw={'height_ratios': height_ratios})
        # fig, axs = plt.subplots(nrows=gm_n + 2, ncols=ncols, figsize=figsize, sharex='col',
        #                        sharey='row', gridspec_kw={'height_ratios': height_ratios})
    else:
        # gm_n = len(genetic_models)
        if figtype == '1':
            figwidth = 8
            figsize = (figwidth * ncols, 5 * gm_n)
        elif figtype == '2':
            figwidth = 5
            figsize = (figwidth * ncols, 2 * gm_n)
        elif figtype == '3':
            # 1.0 for title
            figwidth = 8
            figsize = (figwidth * ncols, 1.0 + 2.5 * gm_n)
        elif figtype == '4':
            figwidth = 8
            figsize = (figwidth * ncols, 4.5 * gm_n)
        elif figtype == '5':
            figwidth = 7
            figsize = (figwidth * ncols, 4 * gm_n)
        elif figtype == '6':
            figwidth = 4.5
            figsize = (figwidth * ncols, 7 * gm_n)
        elif figtype == '9':
            figwidth = 11
            figsize = (figwidth * ncols, 5 * gm_n)
        else:
            raise RuntimeError()

        fig, axs = plt.subplots(nrows=gm_n, ncols=ncols, figsize=figsize, sharex=sharex, sharey=sharey)
        # fig, axs = plt.subplots(nrows=gm_n, ncols=ncols, figsize=figsize, sharex='col', sharey='row')

    return fig, axs, figwidth


# sida to rs
sida_rs = {
    '2:144385104:C:T': 'rs2731561',
    '2:65586769:C:A': 'rs10193337',
    '2:111600519:T:C': 'rs13395354',
    '2:228968345:G:A': 'rs12479220',
    '3:119123814:C:T': 'rs6773050',
    '4:145597759:A:G': 'rs2055059',
    '6:168413393:A:G': 'rs11751451',
    '8:8936749:A:G': 'rs6990464',
    '8:99565446:T:C': 'rs35696573',
    '11:96523529:T:C': 'rs1530117',
    '11:97281245:A:G': 'rs10895739',
    '11:97734972:C:T': 'rs10789728',
    '11:98176280:C:T': 'rs12224719',
    '18:60004253:G:A': 'rs7237982',
    '20:41711226:T:C': 'rs6065574',
    '22:20876509:G:A': 'rs7291930',
    '22:20884910:T:C': 'rs6519408',
}

cb = sns.color_palette("colorblind", 10)
# ['overdom', 'dom','add','rec', 'overrec']
# should be different colors as odds ratio
# mypalette4 = [
#    cb[6], cb[4], cb[2], cb[0], cb[9]
# ]

br = sns.color_palette("bright", 10)
# TMP for reviewer
mypalette4 = [
    # cb[6],
    br[6], cb[4], cb[2], cb[0], cb[9]
]


def plot_manhattan_overlay_pos(
        ss,
        y,
        y_overlay,
        annot=None,  # column name of annotation
        significant=None,
        ymax=None,
        chrom=None,
        ax=None,
        yneg=False,
        gradcolor=None,
        y_name=None,
        overlay_hue=None,
        hue_thre=None,
        annot_thre_lower=None,
        annot_snv_sidas=None,  # either thres_lower or snv_sids
        annot_snv_sidas_rs=None,
        annot_model=None,
        annot_pos=None,  # None, middle
        vmin=None,
        vmax=None,
        xaxis_skip=[],
        xlim=None,
        color=None,
        title=None,
        # x_label_width=None,
        xticks=1000000,
        use_legend=False,
        legend_col=None,
        col_overlay2=None,
        yticks=None,
        yticks_label=None,
        xlabel=None,
):
    """
    # rasterized
    gradcolor is for y_overlay. (if want to use on y, use plot_manhattan_colorbar())
    gradcolor: None, True, (col name): if True, use 'boost_iteration_sprse'

    # x is position

    # TODO: create col='is_boost' to avoid ~ss['boost_alpha'].isna()

    """

    if chrom is None:
        raise NotImplementedError('ny for several chroms')

    if len(pd.unique(ss['chrom'])) > 1:
        raise NotImplementedError('ny for several chroms')

    if (annot_thre_lower is not None) and (annot_snv_sidas is not None):
        raise RuntimeError("cannot indicate both")

    im = None

    if chrom is None:
        if ss[y].isnull().all():
            return ax, im
    else:
        if ss.loc[ss['chrom'] == chrom, y].isnull().all():
            return ax, im

    if color is None:
        color = ['#000088']
    # color_g = ['#646464']
    # color_g = ['#d3d3d3']
    color_g = ['#ebebeb']
    color_hue = ['#880000']

    # print("labels")

    start = int(ss['pos'].iloc[0])
    end = int(ss['pos'].iloc[-1])
    # print('s,e', start, end)

    ss_plot = ss.copy()

    if xlim is True:
        # xlim: (0, last), default add empty before and last
        xmin = ss_plot['pos'].min()
        xmax = ss_plot['pos'].max()
        xlim = (xmin, xmax)

    # print("plot")

    # background
    # plot not in boost_alpha
    if len(ss_plot[ss_plot[y_overlay].isna()]) > 0:
        print("not selected", len(ss_plot[ss_plot[y_overlay].isna()]))
        linewidth = 0
        s = 10
        # s=min(10,1000000/len(ss_plot))
        # s=min(10,10000000/len(ss_plot))

        xs = ss_plot.loc[ss_plot[y_overlay].isna(), 'pos']
        ys = ss_plot.loc[ss_plot[y_overlay].isna(), y]
        ax.scatter(xs, ys, s=s, color=color_g, linewidths=linewidth,
                   label='Not selected', rasterized=True)

        # using pandas interfere ylabel for subplots sharey
        # ax = ss_plot[ss_plot[y_overlay].isna()].plot(kind='scatter', x='pos', y=y,
        #                                             color=color_g, ax=ax, s=s,
        #                                             linewidth=linewidth,
        #                                             rasterized=True,
        #                                             label='Not selected')

    # plot also in boost_alpha
    # plot or not
    if col_overlay2 in ss_plot:
        # plot overlay2 even if boost snvs are not in the region.
        plot_boost = True
    else:
        plot_boost = (len(ss_plot[~ss_plot[y_overlay].isna()]) > 0)
    # if len(ss_plot[~ss_plot[y_overlay].isna()]) > 0:
    if plot_boost:
        print("selected", len(ss_plot[~ss_plot[y_overlay].isna()]))
        s = 200
        # s = 100
        # s_small = 5
        # s=min(20,1000000/len(ss_plot))
        # s=min(20,1000000/len(ss_plot))

        # SNVs in boosting usual plot
        if (gradcolor is None) and (overlay_hue is None):
            ss_plot_overlay = ss_plot[~ss_plot[y_overlay].isna()]
            if ('boost_model_est' in ss_plot_overlay) and (col_overlay2 in ss_plot_overlay):
                color_red = '#E40000'
                color_green = '#00E400'
                markerd = {'overrec': 'P', 'rec': 'X', 'add': 'o', 'dom': '^', 'overdom': 'v'}
                linewidth = 0

                # for finemap
                ss_plot_fine = ss_plot[ss_plot[col_overlay2]]
                print('len fine', len(ss_plot_fine))
                xs = ss_plot_fine['pos']
                ys = ss_plot_fine[y]
                ax.scatter(xs, ys, s=s, color=color_green, linewidths=linewidth,
                           rasterized=True)
                # ax = ss_plot_fine.plot(kind='scatter', x='pos', y=y,
                #                       color=color_green, ax=ax, s=s,
                #                       linewidth=linewidth,
                #                       rasterized=True)
                # for boost
                print('len boost', len(ss_plot_overlay))
                for model in markerd.keys():
                    ss_plot_overlay_model = ss_plot_overlay[ss_plot_overlay['boost_model_est'] == model]
                    # print('len boost', len(ss_plot_overlay_model))
                    xs = ss_plot_overlay_model['pos']
                    ys = ss_plot_overlay_model[y]
                    ax.scatter(xs, ys, s=s, color=color, linewidths=linewidth,
                               marker=markerd[model],
                               rasterized=True)
                    # ax = ss_plot_overlay_model.plot(kind='scatter', x='pos', y=y,
                    #                                color=color, ax=ax, s=s,
                    #                                marker=markerd[model],
                    #                                linewidth=linewidth,
                    #                                rasterized=True)
                # for both
                print('len boost', len(ss_plot_overlay[ss_plot_overlay[col_overlay2]]))
                for model in markerd.keys():
                    ss_plot_overlay_model = ss_plot_overlay[(ss_plot_overlay['boost_model_est'] == model) & (ss_plot_overlay[col_overlay2])]
                    # ss_plot_overlay_model = ss_plot_overlay[(ss_plot_overlay['boost_model_est'] == model) & (ss_plot_overlay['in_finemap'])]
                    # print('len both', len(ss_plot_overlay_model))
                    xs = ss_plot_overlay_model['pos']
                    ys = ss_plot_overlay_model[y]
                    ax.scatter(xs, ys, s=s, color=color_red, linewidths=linewidth,
                               marker=markerd[model],
                               rasterized=True)
                    # ax = ss_plot_overlay_model.plot(kind='scatter', x='pos', y=y,
                    #                                color=color_red, ax=ax, s=s,
                    #                                marker=markerd[model],
                    #                                linewidth=linewidth,
                    #                                rasterized=True)

            elif 'boost_model_est' in ss_plot_overlay:
                # change marker for each model
                markerd = {'add': 'o', 'dom': '^', 'overdom': 'v', 'rec': 'X', 'overrec': 'P', }
                # markerd = {'overrec': 'P', 'rec': 'X', 'add': 'o', 'dom': '^', 'overdom': 'v'}
                markerd_name = {'overrec': 'Overrecessive', 'rec': 'Recessive', 'add': 'Additive', 'dom': 'Dominant', 'overdom': 'Overdominant'}

                markers_for_color = ['overdom', 'dom', 'add', 'rec', 'overrec']
                # markers_for_color = ['add', 'overdom', 'dom', 'overrec', 'rec']
                markerd_color = dict(zip(markers_for_color, mypalette4))
                # markerd = {'overrec': 'X', 'rec': 'x', 'add': '.', 'dom': '+', 'overdom': 'P'}
                # marker edge width
                linewidth = 0
                # for model in markerd.keys():
                for model in markers_for_color:
                    ss_plot_overlay_model = ss_plot_overlay[ss_plot_overlay['boost_model_est'] == model]

                    # print('ss_plot_overlay_model model', model, ss_plot_overlay_model)
                    # print('p', model, ss_plot_overlay_model.loc[:, '-log10(P)'])
                    xs = ss_plot_overlay_model['pos']
                    ys = ss_plot_overlay_model[y]
                    ax.scatter(xs, ys, s=s,
                               color=markerd_color[model],
                               linewidths=linewidth,
                               marker=markerd[model],
                               label=markerd_name[model],
                               rasterized=True)

                    # ax = ss_plot_overlay_model.plot(kind='scatter', x='pos', y=y,
                    #                                color=markerd_color[model],
                    #                                # color=color,
                    #                                ax=ax, s=s,
                    #                                marker=markerd[model],
                    #                                linewidth=linewidth,
                    #                                rasterized=True,
                    #                                label=markerd_name[model])

            else:
                # ax = ss_plot[~ss_plot[y_overlay].isna()].plot(kind='scatter', x='serial', y=y,
                xs = ss_plot_overlay['serial']
                ys = ss_plot_overlay[y]
                ax.scatter(xs, ys, s=s,
                           color=color,
                           rasterized=True)
                # ax = ss_plot_overlay.plot(kind='scatter', x='serial', y=y,
                #                          color=color, ax=ax, s=s,
                #                          rasterized=True)

            if annot is not None:
                # annot indicated indexs only and change plot color
                if annot_snv_sidas is not None:
                    def plot_text(ax, ss_plot_overlay, i, y):
                        x_pos = ss_plot_overlay['pos'].iloc[i]
                        y_pos = ss_plot_overlay[y].iloc[i]
                        y_annot = 20
                        # y_annot=y_pos+20+(i%3)*5
                        # y_annot = 20 + (j % 4) * 5
                        # avoid overlap
                        if ss_plot_overlay[annot].iloc[i] == '11:120222475:T:C':
                            x_pos_txt = x_pos + 500000
                        else:
                            x_pos_txt = x_pos

                        ax.text(x_pos_txt, y_annot,
                                ss_plot_overlay[annot].iloc[i],
                                color='black',
                                horizontalalignment='center',
                                fontsize='x-small'
                                )
                        # draw line
                        ax.plot([x_pos, x_pos], [y_pos + 1, y_annot - 1], linewidth=1,
                                color='#696969')
                        # ax.plot(x_pos,y_pos,x_pos,y_annot)
                        # ax.plot(ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y].iloc[i],ss_plot_overlay['serial'].iloc[i],y_annot,)
                        # linewidth
                        print("i", i)
                        # print("a",ss_plot_overlay['serial'])
                        # print("serial",ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y]+10,
                        # ss_plot_overlay['sida'].iloc[i])
                        return ax

                    # palette=sns.color_palette("hls", 8)
                    # palette=sns.color_palette("husl", 8)
                    # for palette
                    color_a = '#E40000'
                    linewidth = 0
                    for i in range(ss_plot_overlay.shape[0]):
                        if ss_plot_overlay['sida'].iloc[i] in annot_snv_sidas:
                            print("annot", i, ss_plot_overlay['sida'].iloc[i])
                            ax = plot_text(ax, ss_plot_overlay, i, y)

                            xs = ss_plot_overlay.iloc[[i]]['pos']
                            ys = ss_plot_overlay.iloc[[i]][y]
                            ax.scatter(xs, ys, s=s,
                                       color=color_a,
                                       linewidths=linewidth,
                                       marker=markerd[ss_plot_overlay['boost_model_est'].iloc[i]],
                                       rasterized=True)

                            # ax = ss_plot[~ss_plot[y_overlay].isna()].plot(kind='scatter', x='serial', y=y,
                            # ax = ss_plot_overlay.iloc[[i]].plot(kind='scatter', x='pos', y=y,
                            #                                    color=color_a, ax=ax, s=s,
                            #                                    marker=markerd[ss_plot_overlay['boost_model_est'].iloc[i]],
                            #                                    linewidth=linewidth,
                            #                                    rasterized=True)

                # annot indicated by 'sida' and write 'rs' but do not change plot color
                elif annot_snv_sidas_rs is not None:
                    print('annot_snv_sidas_rs', annot)

                    def plot_text(ax, ss_plot_overlay, i, y, annot_pos):

                        varname = sida_rs[ss_plot_overlay[annot].iloc[i]]
                        # varname= ss_plot_overlay[annot].iloc[i]

                        x_pos = ss_plot_overlay['pos'].iloc[i]
                        y_pos = ss_plot_overlay[y].iloc[i]

                        if annot_pos is None:
                            # FIXME
                            # ideally, want to use y_pos+ylim*0.05
                            # but ylim might change later due to sharey
                            # ywid = ss_plot_overlay[y].max() - ss_plot_overlay[y].min()
                            ywid = ax.get_ylim()[1] - ax.get_ylim()[0]
                            if ss_plot_overlay[annot].iloc[i] == '2:65586769:C:A':
                                y_annot = y_pos + ywid * 0.15
                            elif y == 'decreasing_loss' and ss_plot_overlay[annot].iloc[i] == '11:98176280:C:T':
                                y_annot = y_pos - ywid * 0.05
                            else:
                                y_annot = y_pos + ywid * 0.05
                            # y_annot = y_pos*1.1
                            # y_annot = y_pos
                            # y_annot = 20
                            # y_annot=y_pos+20+(i%3)*5
                            # y_annot = 20 + (j % 4) * 5
                            # avoid overlap
                            # if ss_plot_overlay[annot].iloc[i] == '11:120222475:T:C':
                            #    x_pos_txt = x_pos #+ 500000
                            # else:
                            x_pos_txt = x_pos + 50000

                            # bbox={'facecolor':'white',
                            #      'alpha':0.5}
                            ax.text(x_pos_txt, y_annot,
                                    varname,
                                    color='black',
                                    horizontalalignment='left',
                                    # horizontalalignment='center',
                                    verticalalignment='center',
                                    fontsize='x-small',
                                    # bbox=bbox
                                    )
                            # draw line
                            ax.plot([x_pos, x_pos_txt], [y_pos, y_annot], linewidth=1,
                                    color='#696969',
                                    zorder=0,
                                    )
                        elif annot_pos == 'middle':
                            ywid = ax.get_ylim()[1] - ax.get_ylim()[0]
                            x_pos_txt = x_pos + 50000

                            if ss_plot_overlay[annot].iloc[i] == '2:65586769:C:A':
                                y_annot = y_pos + ywid * 0.5
                                ax.text(x_pos_txt, y_annot,
                                        varname,
                                        color='black',
                                        horizontalalignment='center',
                                        verticalalignment='bottom',
                                        fontsize='x-small',
                                        # bbox=bbox
                                        )
                            elif y == 'decreasing_loss' and ss_plot_overlay[annot].iloc[i] == '11:98176280:C:T':
                                y_annot = y_pos - ywid * 0.02
                                ax.text(x_pos_txt, y_annot,
                                        varname,
                                        color='black',
                                        horizontalalignment='center',
                                        verticalalignment='top',
                                        fontsize='x-small',
                                        # bbox=bbox
                                        )
                            else:
                                y_annot = y_pos + ywid * 0.02

                                ax.text(x_pos_txt, y_annot,
                                        varname,
                                        color='black',
                                        horizontalalignment='center',
                                        verticalalignment='bottom',
                                        fontsize='x-small',
                                        # bbox=bbox
                                        )
                            # draw line
                            ax.plot([x_pos, x_pos_txt], [y_pos, y_annot], linewidth=1,
                                    color='#696969',
                                    zorder=0,
                                    )
                        else:
                            raise NotImplementedError

                        # ax.plot(x_pos,y_pos,x_pos,y_annot)
                        # ax.plot(ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y].iloc[i],ss_plot_overlay['serial'].iloc[i],y_annot,)
                        # linewidth
                        # print("i", i)
                        # print("a",ss_plot_overlay['serial'])
                        # print("serial",ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y]+10,
                        # ss_plot_overlay['sida'].iloc[i])
                        return ax

                    # palette=sns.color_palette("hls", 8)
                    # palette=sns.color_palette("husl", 8)
                    # for palette
                    color_a = '#E40000'
                    linewidth = 0
                    for i in range(ss_plot_overlay.shape[0]):
                        if ss_plot_overlay['sida'].iloc[i] in annot_snv_sidas_rs:
                            print("annot", i, ss_plot_overlay['sida'].iloc[i])
                            ax = plot_text(ax, ss_plot_overlay, i, y, annot_pos)

                            # ax = ss_plot_overlay.iloc[[i]].plot(kind='scatter', x='pos', y=y,
                            #                                    color=color_a, ax=ax, s=s,
                            #                                    marker=markerd[ss_plot_overlay['boost_model_est'].iloc[i]],
                            #                                    linewidth=linewidth,
                            #                                    rasterized=True)

                # annot all under annot_thre_lower
                elif annot_thre_lower is not None:
                    def plot_text(ax, ss_plot_overlay, i, y, j):
                        x_pos = ss_plot_overlay['pos'].iloc[i]
                        y_pos = ss_plot_overlay[y].iloc[i]
                        # y_annot=y_pos+20+(i%3)*5
                        y_annot = 20 + (j % 4) * 5
                        ax.text(x_pos, y_annot,
                                ss_plot_overlay[annot].iloc[i],
                                color='black',
                                horizontalalignment='center',
                                fontsize='x-small'
                                )
                        # draw line
                        ax.plot([x_pos, x_pos], [y_pos + 1, y_annot], linewidth=1,
                                color='#696969')
                        # ax.plot(x_pos,y_pos,x_pos,y_annot)
                        # ax.plot(ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y].iloc[i],ss_plot_overlay['serial'].iloc[i],y_annot,)
                        # linewidth
                        print("i", i)
                        # print("a",ss_plot_overlay['serial'])
                        # print("serial",ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y]+10,
                        # ss_plot_overlay['sida'].iloc[i])
                        return ax

                    # count for annotation pos
                    j = 0
                    for i in range(ss_plot_overlay.shape[0]):
                        if annot_thre_lower is not None and ss_plot_overlay[y].iloc[i] > annot_thre_lower:
                            continue
                        # annot_model should be None or iteratable
                        if annot_model is not None and ss_plot_overlay['boost_model_est'].iloc[i] not in annot_model:
                            continue
                        ax = plot_text(ax, ss_plot_overlay, i, y, j)
                        j += 1

                # faster
                # if annot_thre_lower is None:
                #    for i in range(ss_plot_overlay.shape[0]):
                #        ax=plot_text(ax,ss_plot_overlay,i,y)
                #        #ax.text(ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y]+10,
                #        #ss_plot_overlay['sida'].iloc[i],
                #        #color='black')
                # else:
                #    ss_plot_overlay_thre=ss_plot_overlay[ss_plot_overlay[y]<annot_thre_lower]
                #    for i in range(ss_plot_overlay_thre.shape[0]):
                #        ax=plot_text(ax,ss_plot_overlay_thre,i,y)
                #        #if ss_plot_overlay[y].iloc[i]>annot_thre_lower:
                #        #    continue
                #        #ax.text(ss_plot_overlay['serial'].iloc[i],ss_plot_overlay[y]+10,
                #        #ss_plot_overlay['sida'].iloc[i],
                #        #color='black',
                #        #horizontalalignment='center',
                #        #)

        if overlay_hue is not None:
            df_plot = ss_plot[~ss_plot[y_overlay].isna()]

            xs = df_plot[df_plot[overlay_hue] >= hue_thre]['pos']
            ys = df_plot[df_plot[overlay_hue] >= hue_thre][y]
            ax.scatter(xs, ys, s=s,
                       color=color_hue,
                       label=overlay_hue + '>=' + str(hue_thre),
                       rasterized=True)
            xs = df_plot[df_plot[overlay_hue] < hue_thre]['pos']
            ys = df_plot[df_plot[overlay_hue] < hue_thre][y]
            ax.scatter(xs, ys, s=s,
                       color=color_hue,
                       label=overlay_hue + '<' + str(hue_thre),
                       rasterized=True)

            # ax = df_plot[df_plot[overlay_hue] >= hue_thre].plot(kind='scatter', x='pos', y=y,
            #                                                    color=color_hue, ax=ax, s=s,
            #                                                    label=overlay_hue + '>=' + str(hue_thre),
            #                                                    rasterized=True
            #                                                    )
            # ax = df_plot[df_plot[overlay_hue] < hue_thre].plot(kind='scatter', x='pos', y=y,
            #                                                   color=color, ax=ax, s=s,
            #                                                   label=overlay_hue + '<' + str(hue_thre),
            #                                                   rasterized=True
            #                                                   )
            # TODO: add text in plot here: https://www.python-graph-gallery.com/46-add-text-annotation-on-scatterplot

        if gradcolor is not None:
            if gradcolor is True:
                gradcolor_ = 'boost_iteration_sparse'
                if gradcolor_ not in ss_plot:
                    ss_plot.loc[:, gradcolor_] = (
                        ss_plot['boost_iteration'] / 100).astype(np.int) * 100
            else:
                # ex. color = 'boost_iteration'
                gradcolor_ = gradcolor

            df_plot = ss_plot[~ss_plot[y_overlay].isna()]
            im = ax.scatter(df_plot['pos'], df_plot[y],
                            c=df_plot[gradcolor_], cmap='jet', s=s,
                            vmin=vmin, vmax=vmax,
                            rasterized=True
                            )

    x_labels = []
    x_labels_pos = []
    # x_labels = [start, end]
    # x_labels_pos = [start, end]
    # should be %1M=0
    # -> allow 500kbp
    # if x_label_width is None:
    #    x_label_width = 2000000
    if xticks is not None:
        x_label_width = xticks
        start_label = ((start + x_label_width - 1) // x_label_width) * x_label_width
        end_label = (end // x_label_width) * x_label_width
        print('s,e', start_label, end_label)
        for label_pos in range(start_label, end_label + x_label_width, x_label_width):
            x_labels_pos.append(label_pos)
            if label_pos % 1000000 == 0:
                # make integer
                x_labels.append(str(label_pos // 1000000) + 'M')
            else:
                x_labels.append(str(label_pos / 1000000) + 'M')

        # print("x_labels_pos", x_labels_pos)
        # print("x_labels", x_labels)
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)

    if xticks == 1000000:
        # set minor ticker
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    elif xticks == 500000:
        # set minor ticker
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    if chrom is None:
        # temporary
        # ax.set_xlim([0, len(ss)])
        ax.set_xlim((0, ss['serial'].max()))
    else:
        # ??
        # ax.set_xlim([0, len(group)])
        pass

    if ymax == '0_default':
        ax.set_ylim(bottom=0)
    else:
        if ymax is None:
            if chrom is None:
                ymax = ss[y].max(skipna=True) * 1.1  # same as default
                # ymax = ss[y].max(skipna=True) * 1.05
                # ymax = ss[y].max(numeric_only=True)*1.05
            else:
                ymax = ss.loc[ss['chrom'] == chrom, y].max(skipna=True) * 1.1
                # ymax = ss.loc[ss['chrom'] == chrom, y].max(skipna=True) * 1.05
                # ymax = ss.loc[ss['chrom'] == chrom, y].max()*1.05
        else:
            ymax = ymax
            print("ymax vs max", ymax, ss[y].max(skipna=True))

        if not yneg:
            ax.set_ylim([0, ymax])
        else:
            ax.set_ylim([-ymax, ymax])

    if yticks is not None:
        if yticks_label is not None:
            kwargs = dict(labels=yticks_label)
        else:
            kwargs = {}
        ax.set_yticks(yticks, minor=False, **kwargs)

    if xlim is not None:
        ax.set_xlim(xlim)
    # TEMP not to indicate fontsize to learn default size
    if xlabel is None:
        ax.set_xlabel(None)
        # ax.set_xlabel('')
    elif xlabel is True:
        ax.set_xlabel('Position in chromosome ' + str(chrom) + ' [bp]')
    else:
        ax.set_xlabel(xlabel)

    if y_name is None:
        ax.set_ylabel(y)
    else:
        ax.set_ylabel(y_name)

    if title is not None:
        ax.set_title(title)

    if significant is not None:
        ax.axhline(y=significant, color='gray', linestyle='-', linewidth=2.0)

    if (overlay_hue is not None) or use_legend:

        # ax.legend()
        # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # TODO: refine
        # for genetic models change order
        handles, labels = plt.gca().get_legend_handles_labels()
        if len(labels) == 6:
            if legend_col is None:
                order = [1, 2, 3, 4, 5, 0]
                # plt.legend([handles[i] for i in order], [labels[i] for i in order])
                ax.legend([handles[i] for i in order], [labels[i] for i in order], bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                print('legend_col=2')
                order = [1, 2, 3, 4, 5, 0]
                # plt.legend([handles[i] for i in order], [labels[i] for i in order])
                # ncol no effects here
                ax.legend([handles[i] for i in order], [labels[i] for i in order], bbox_to_anchor=(1.05, 1), loc='upper left', ncol=legend_col)

        else:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        if legend_col is None:
            pcommon.move_legend_ax(ax, title='Genetic model')
        else:
            pcommon.move_legend_ax(ax, title='Genetic model', ncol=legend_col)

    else:
        pass
        # if not using pandas, legend will not be generated automatically
        # ax.get_legend().remove()

    # ax.set_xticks(fontsize=32, rotation=60)
    # ax.set_yticks(fontsize=32)
    # plt.axhline(y=significant, color='gray', linestyle='-', linewidth = 2.0)
    # plt.xticks(fontsize=32, rotation=60)
    # plt.yticks(fontsize=32)

    # return fig,ax
    return ax, im


def manhattan_loss_multi(ssobs_wgt_, ax, title, annot, chrom, color, x_label_width, ymax, annot_pos=None):
    y_name = 'Non-additive GenoBoost\'s\npredictive utility'
    # y_name = 'Decreasing loss'

    yticks = None
    yticks_label = None
    # for t2d
    # yticks = [0.0, 0.0001, 0.0002, 0.0003]
    # yticks_label = ['0.0', '1e-4', '2e-4', '3e-4']
    # ylim=[0,ymax] since yneg=False
    ax = plot_manhattan_overlay_pos(ssobs_wgt_, y='decreasing_loss',
                                    y_overlay='boosting_nonadd_effect_abs',
                                    # y_overlay='boost_alpha',
                                    ax=ax, xlim=True,
                                    chrom=chrom, color=color, title=title,
                                    xticks=x_label_width,
                                    y_name=y_name,
                                    yticks=yticks,
                                    yticks_label=yticks_label,
                                    annot='sida',
                                    # annot='boosting_var',
                                    annot_snv_sidas_rs=annot,
                                    annot_pos=annot_pos,
                                    ymax=ymax,
                                    )
    return ax


def manhattan_man_multi(ssobs_wgt_, ax, genetic_model, title, annot, chrom, color, x_label_width, ymax, xlabel, annot_pos=None):
    if (genetic_model == 'add') or (genetic_model == 'condadd'):
        y_name = '-log($P_{add}$)'
    elif (genetic_model == 'dom') or (genetic_model == 'conddom'):
        y_name = '-log($P_{dom}$)'
    elif (genetic_model == 'rec') or (genetic_model == 'condrec'):
        y_name = '-log($P_{rec}$)'
    elif (genetic_model == 'hetonly') or (genetic_model == 'condhetonly'):
        y_name = '-log($P_{hetonly}$)'
    else:
        y_name = None

    # print("p to plot", ssobs_wgt_['-log10(P)'])

    ax = plot_manhattan_overlay_pos(ssobs_wgt_, y='-log10(P)',
                                    y_overlay='boosting_nonadd_effect_abs',
                                    # y_overlay='boost_alpha',
                                    ax=ax, xlim=True,
                                    chrom=chrom, color=color, title=title,
                                    xticks=x_label_width,
                                    y_name=y_name,
                                    annot='sida',
                                    # TMP
                                    # annot='boost_var',
                                    annot_snv_sidas_rs=annot,
                                    annot_pos=annot_pos,
                                    ymax=ymax,
                                    xlabel=xlabel
                                    )
    return ax


# TMP
def add_sidv(df):
    # df['sidv'] = df['chrom'].astype(
    df.loc[:, 'sidv'] = df['chrom'].astype(
        str) + ':' + df['pos'].astype(str) + ':' + df['a2'] + '>' + df['a1']
    return df


def get_genetic_modeld():
    return {'overrec': 'overrecessive', 'rec': 'recessive', 'add': 'additive', 'dom': 'dominant', 'overdom': 'overdominant'}


def get_rowi(start_text, rowd_end):
    padding_width = 0

    for rowi in range(len(rowd_end)):
        row_end = rowd_end[rowi]
        # print(row_end , padding_width , start_text)

        if row_end + padding_width < start_text:
            return rowi
    rowi = len(rowd_end)
    return rowi


def text_pos(gene_name, start, width_x_per_letter):
    # for ax.text()
    # padding for 1/4 letter
    text_end = start - width_x_per_letter * 0.25

    text_width = len(gene_name) * width_x_per_letter
    # print('text_end', text_end, text_width)
    # just for rough estimation
    text_start = text_end - text_width
    return text_end, text_start


def compute_row_max(fontsize, width_inch, ax, gene_exon):
    # width_ax_pt = width_inch * 72
    # xlims = ax.get_xlim()
    # width_ax_x = xlims[1] - xlims[0]
    # width_x_per_pt = width_ax_x / width_ax_pt
    # width_x_per_letter = width_x_per_pt * fontsize

    width_x_per_letter = get_width_x_per_letter(width_inch, ax, fontsize)

    # TODO: rename
    (xlim_lower, xlim_upper) = ax.get_xlim()

    # get row_max
    rowd_end = {}  # dict of limit of x position
    for gene_id in gene_exon['gene_id'].unique():
        # gene
        gene_exon_id = gene_exon[gene_exon['gene_id'] == gene_id]
        # assume gene_id is unique
        gene_ = gene_exon_id[gene_exon_id['gene_type'] == 'gene'].iloc[0]
        gene_name = gene_.loc['gene']
        start = max(gene_.loc['start'], xlim_lower)
        end = min(gene_.loc['end'], xlim_upper)
        # start =gene_.loc['start'],xlim_begin
        # end = gene_.loc['end']
        # width = end - start

        # text_width=(len(gene_name)+1)*width_x_per_letter
        # start_text=start-text_width
        _, text_start = text_pos(gene_name, start, width_x_per_letter)

        rowi = get_rowi(text_start, rowd_end)
        rowd_end[rowi] = end

    row_max = len(rowd_end)
    # print('row_max', row_max)
    return row_max


def format_gene_exon(gene_exon, chrompos, ax):
    (xlim_lower, xlim_upper) = ax.get_xlim()
    chrom, _, _ = chrompos
    filt = (gene_exon['chrom'] == chrom) & (gene_exon['end'] >= xlim_lower) & (gene_exon['start'] <= xlim_upper)

    # usually xlim is narrower than chrompos and some genes are in chrompos but not in xlim
    # chrom, pos_start, pos_end = chrompos
    # if start or end is in the region
    # filt = (gene_exon['chrom'] == chrom) & (gene_exon['end'] >= pos_start) & (gene_exon['start'] <= pos_end)
    # filt = (gene_exon['chrom'] == chrom) & (gene_exon['start'] >= pos_start) & (gene_exon['end'] <= pos_end)
    gene_exon_chrompos = gene_exon.loc[filt, :].copy()
    return gene_exon_chrompos


def search_best_fontsize(fontsize_start, width_inch, row_max_thre, ax, gene_exon):
    # adjust fontsize so taht row_max<row_thre

    fontsize = fontsize_start
    while True:
        # print('search next fontsize: ', fontsize)
        row_max = compute_row_max(fontsize, width_inch, ax, gene_exon)
        if row_max <= row_max_thre:
            return fontsize, row_max
        fontsize -= 1.0

        if fontsize < 3.0:
            raise RuntimeError('Too small fontsize: give me more space')


def get_width_x_per_letter(width_inch, ax, fontsize):
    width_ax_pt = width_inch * 72
    xlims = ax.get_xlim()
    width_ax_x = xlims[1] - xlims[0]
    width_x_per_pt = width_ax_x / width_ax_pt
    width_x_per_letter = width_x_per_pt * fontsize
    return width_x_per_letter


def plot_exon(ax, gene_exon, width_inch):
    # how to get width_inch from ax??
    # print('ax size',ax.get_position().bounds)
    # l, b, w, h = ax.get_position().bounds

    # pt
    fontsize_start = 17
    row_max_thre = 8

    # print('xlim', ax.get_xlim())

    fontsize, row_max = search_best_fontsize(fontsize_start, width_inch, row_max_thre, ax, gene_exon)
    # print('use fontsize', fontsize)

    width_x_per_letter = get_width_x_per_letter(width_inch, ax, fontsize)

    # get row_max
    # rowd_end = {}  # dict of limit of x position
    # for gene_id in gene_exon['gene_id'].unique():
    #    # gene
    #    gene_exon_id = gene_exon[gene_exon['gene_id'] == gene_id]
    #    # assume gene_id is unique
    #    gene_ = gene_exon_id[gene_exon_id['gene_type'] == 'gene'].iloc[0]
    #    gene_name = gene_.loc['gene']
    #    start = gene_.loc['start']
    #    end = gene_.loc['end']
    #    width = end - start

    #    # text_width=(len(gene_name)+1)*width_x_per_letter
    #    # start_text=start-text_width
    #    _, text_start = text_pos(gene_name, start, width_x_per_letter)

    #    rowi = get_rowi(text_start, rowd_end)
    #    rowd_end[rowi] = end

    # row_max = len(rowd_end)
    # print('row_max', row_max)

    ax.set_ylim((0.0, 1.0))
    # ylim is 0-1.0
    row_height = min((1.0 / row_max), 0.2)
    # row_height=min((0.9/row_max),0.2)
    height_exon = row_height * 0.5
    height_gene = row_height * 0.1

    (xlim_lower, xlim_upper) = ax.get_xlim()
    print('xlim', ax.get_xlim())

    # rowi=0
    rowd_end = {}  # dict of limit of x position
    for gene_id in gene_exon['gene_id'].unique():
        # print('gene', gene_id)

        # gene
        gene_exon_id = gene_exon[gene_exon['gene_id'] == gene_id]
        # print('gene_exon_id', gene_exon_id)
        # assume gene_id is unique
        gene_ = gene_exon_id[gene_exon_id['gene_type'] == 'gene'].iloc[0]
        gene_name = gene_.loc['gene']
        start = max(gene_.loc['start'], xlim_lower)
        end = min(gene_.loc['end'], xlim_upper)
        # start = gene_.loc['start']
        # end = gene_.loc['end']
        width = end - start

        # +1 is for padding
        text_width = (len(gene_name) + 1) * width_x_per_letter
        text_start = start - text_width

        text_end, text_start = text_pos(gene_name, start, width_x_per_letter)

        rowi = get_rowi(text_start, rowd_end)

        black = '#000000'
        # 0.5 for center
        ycenter = 1.0 - (rowi + 0.5) * row_height
        # ycenter=0.9 - (rowi+0.5)*row_height
        # ytop=0.9 - rowi*row_height
        # height=0.2
        lw = 1.0
        r = patches.Rectangle(xy=(start, ycenter - height_gene / 2), width=width, height=height_gene, fc=black, ec=black, fill=True, lw=lw)
        # r = patches.Rectangle(xy=(start, 0.8), width=width, height=-0.5, fc=black,ec=black, fill=True)
        ax.add_patch(r)

        for _, exon_row in gene_exon_id[gene_exon_id['gene_type'] == 'exon'].iterrows():
            # do not use max( , xlim_lower)
            start = exon_row.loc['start']
            end = exon_row.loc['end']
            width = end - start
            # ng: linewidth=0.0
            r = patches.Rectangle(xy=(start, ycenter - height_exon / 2), width=width, height=height_exon, fc=black, ec=black, fill=True, lw=lw)
            ax.add_patch(r)

        # ycenter=ytop-height/2
        ax.text(x=text_end, y=ycenter, s=gene_name, va='center', ha='right', fontsize=fontsize)
        # ax.text(x=start_text,y=ytop,s=gene_name,va='top',ha='right',fontsize=fontsize)

        rowd_end[rowi] = end

        # exon

    return ax


def plot_manhattan_multi(ssd,
                         gene_exon,
                         phe_title,
                         chrompos,
                         genetic_models,
                         legend, plot_loss, plot_gene,
                         fig, axs, figwidth, annot=None, figtype='1', ymax=None):

    print('axs', axs.shape)

    chrom = chrompos[0]
    # chrom, pos_start, pos_end = chrompos

    def format_ssobs(ssobs_wgt, chrompos):
        # TMP comment
        # create 'boost_iteration_sparse'
        # ssobs_wgt.loc[:, 'boost_iteration_sparse'] = (ssobs_wgt['boost_iteration'] / 100).astype(np.int) * 100

        ssobs_wgt = add_sidv(ssobs_wgt)
        if 'boost_model_est' not in ssobs_wgt:
            raise RuntimeError('TMP')
            ssobs_wgt.loc[:, 'boost_model_est'] = 'add'
        ssobs_wgt.loc[:, 'Boost_model_est'] = ssobs_wgt['boost_model_est'].map(get_genetic_modeld())
        ssobs_wgt.loc[:, 'annot'] = ssobs_wgt['sidv'].astype(str) + '\n(' + ssobs_wgt['Boost_model_est'].astype(str) + ')'

        # extract chrompos
        chrom, pos_start, pos_end = chrompos
        filt = (ssobs_wgt['chrom'] == chrom) & (ssobs_wgt['pos'] >= pos_start) & (ssobs_wgt['pos'] <= pos_end)
        ssobs_wgt = ssobs_wgt.loc[filt, :].copy()
        # print('chrompos len', len(ssobs_wgt))

        if len(ssobs_wgt) == 0:
            raise RuntimeError('len=0')
        return ssobs_wgt

    color = ['#dc0000']

    ssobs_wgt_ = ssd['add']
    # logger.debug('ssobs_wgt ori {}'.format(ssobs_wgt_))
    logger.debug('ssobs_wgt ori {}'.format(ssobs_wgt_.columns))
    ssobs_wgt_ = format_ssobs(ssobs_wgt_, chrompos)
    xmin = ssobs_wgt_['pos'].min()
    xmax = ssobs_wgt_['pos'].max()
    if xmax - xmin > 2000000:
        x_label_width = 1000000
    else:
        x_label_width = 500000

    # x_label_width = 1000000
    # x_label_width = 500000
    # if (phe == 'bc') or (phe == 'af'):
    #    x_label_width = 1000000
    # else:
    #    x_label_width = 2000000

    if plot_loss:
        genetic_models = ('loss',) + genetic_models
    # print('genetic models', genetic_models)

    if plot_gene:
        axs_man = axs[:-2]
    else:
        axs_man = axs

    assert len(axs_man) == len(genetic_models)

    for i, gmodel in enumerate(genetic_models):
        # if len(genetic_models)==1:
        #    ax=axs
        # else:
        #    ax=axs[i]
        ax = axs[i]

        # or set title for all ax and later hide them
        if i == 0:
            # print('annot', annot)
            # assert len(annot) == 1
            if (annot[0] is None) or len(annot) > 1:
                # title = validate_paper.get_phed()[phe]
                title = phe_title
            else:
                sida = annot[0]
                # TMP
                rs = sida_rs[sida]
                title = phe_title + ': ' + rs
            # title = validate_paper.get_phed()[phe] + '\nchromosome ' + str(chrom)
            # title = validate_paper.get_phed()[phe] + ' chromosome ' + str(chrom)
        else:
            title = ''

        if gmodel == 'loss':
            # if 'loss' not in ssobs_wgt_.columns:
            #    print('skip since loss does not exist.')
            #    return None

            # loss is in any ssobs
            ssobs_wgt_ = ssd['add'].copy()
            # print('loss in?', ssobs_wgt_.columns)
            ssobs_wgt_ = format_ssobs(ssobs_wgt_, chrompos)
            manhattan_loss_multi(ssobs_wgt_, ax, title, annot, chrom, color, x_label_width, ymax)

        else:
            # print('ssd', ssd.keys())
            ssobs_wgt_ = ssd[gmodel].copy()
            ssobs_wgt_ = format_ssobs(ssobs_wgt_, chrompos)
            manhattan_man_multi(ssobs_wgt_, ax, gmodel, title, annot, chrom, color, x_label_width, ymax, xlabel=True)

    # show x label for last man plot
    axs_man[-1].tick_params(axis='x', labelbottom=True)
    # axs_man[-1].tick_params(axis='x', labelbottom=True)
    # ax.tick_params(axis='x', labelbottom=True)

    # set xlabel only for the last man plot
    for ax in axs_man[:-1]:
        ax.xaxis.set_tick_params(labelbottom=False)
        # ax.tick_params(axis='x',labelbottom=False)
        ax.set_xlabel(None)

    if plot_gene:
        # plot nothing on axs[-2]
        axs[-2].axis('off')

        # TMP
        ax = axs[-1]
        gene_exon_ = format_gene_exon(gene_exon, chrompos, ax)
        print('xlim', ax.get_xlim())
        # plot gene
        plot_exon(ax, gene_exon_, figwidth)
        # plot_exon(ax, gene_exon_, figsize[0])
        # plot_exon(axs[-1],gene_exon_)
        # ax.tick_params(axis='both',labelbottom=False,labelleft=False)

        ax.axis('off')

    return fig


# mainly for slide
def plot_manhattan_cols(ssd,
                        gene_exon,
                        phe_title,
                        chrompos,
                        genetic_models,
                        legend, plot_loss, plot_gene,
                        fig, axs, figwidth, annot=None, figtype='1', ymax=None,
                        annot_pos=None):

    # raise NotImplementedError

    print('axs', axs.shape)

    chrom = chrompos[0]
    # chrom, pos_start, pos_end = chrompos

    def format_ssobs(ssobs_wgt, chrompos):
        # TMP comment
        # create 'boost_iteration_sparse'
        # ssobs_wgt.loc[:, 'boost_iteration_sparse'] = (ssobs_wgt['boost_iteration'] / 100).astype(np.int) * 100

        ssobs_wgt = add_sidv(ssobs_wgt)
        if 'boost_model_est' not in ssobs_wgt:
            raise RuntimeError('TMP')
            ssobs_wgt.loc[:, 'boost_model_est'] = 'add'
        ssobs_wgt.loc[:, 'Boost_model_est'] = ssobs_wgt['boost_model_est'].map(get_genetic_modeld())
        ssobs_wgt.loc[:, 'annot'] = ssobs_wgt['sidv'].astype(str) + '\n(' + ssobs_wgt['Boost_model_est'].astype(str) + ')'

        # extract chrompos
        chrom, pos_start, pos_end = chrompos
        filt = (ssobs_wgt['chrom'] == chrom) & (ssobs_wgt['pos'] >= pos_start) & (ssobs_wgt['pos'] <= pos_end)
        ssobs_wgt = ssobs_wgt.loc[filt, :].copy()
        # print('chrompos len', len(ssobs_wgt))

        if len(ssobs_wgt) == 0:
            raise RuntimeError('len=0')
        return ssobs_wgt

    color = ['#dc0000']

    ssobs_wgt_ = ssd['add']
    # logger.debug('ssobs_wgt ori {}'.format(ssobs_wgt_))
    logger.debug('ssobs_wgt ori {}'.format(ssobs_wgt_.columns))
    ssobs_wgt_ = format_ssobs(ssobs_wgt_, chrompos)
    xmin = ssobs_wgt_['pos'].min()
    xmax = ssobs_wgt_['pos'].max()
    if xmax - xmin > 2000000:
        x_label_width = 1000000
    else:
        x_label_width = 500000

    # x_label_width = 1000000
    # x_label_width = 500000
    # if (phe == 'bc') or (phe == 'af'):
    #    x_label_width = 1000000
    # else:
    #    x_label_width = 2000000

    if plot_loss:
        genetic_models = ('loss',) + genetic_models
    # print('genetic models', genetic_models)

    if plot_gene:
        raise NotImplementedError
        axs_man = axs[:-2]
    else:
        axs_man = axs

    assert len(axs_man) == len(genetic_models)

    for i, gmodel in enumerate(genetic_models):
        # if len(genetic_models)==1:
        #    ax=axs
        # else:
        #    ax=axs[i]
        ax = axs[i]

        # FIXME: show 'GWAS(additive)
        # if i == 0:
        #    # print('annot', annot)
        #    # assert len(annot) == 1
        #    if (annot[0] is None) or len(annot) > 1:
        #        # title = validate_paper.get_phed()[phe]
        #        title = phe_title
        #    else:
        #        sida = annot[0]
        #        # TMP
        #        rs = sida_rs[sida]
        #        title = phe_title + ': ' + rs
        #    # title = validate_paper.get_phed()[phe] + '\nchromosome ' + str(chrom)
        #    # title = validate_paper.get_phed()[phe] + ' chromosome ' + str(chrom)
        # else:
        # title = ''
        title = ''

        if gmodel == 'loss':
            # if 'loss' not in ssobs_wgt_.columns:
            #    print('skip since loss does not exist.')
            #    return None

            # loss is in any ssobs
            ssobs_wgt_ = ssd['add'].copy()
            # print('loss in?', ssobs_wgt_.columns)
            ssobs_wgt_ = format_ssobs(ssobs_wgt_, chrompos)
            manhattan_loss_multi(ssobs_wgt_, ax, title, annot, chrom, color, x_label_width, ymax, annot_pos=annot_pos)

        else:
            # print('ssd', ssd.keys())
            ssobs_wgt_ = ssd[gmodel].copy()
            ssobs_wgt_ = format_ssobs(ssobs_wgt_, chrompos)
            manhattan_man_multi(ssobs_wgt_, ax, gmodel, title, annot, chrom, color, x_label_width, ymax,
                                xlabel=None, annot_pos=annot_pos)
            # xlabel=True)

    for ax in axs_man:
        ax.xaxis.set_tick_params(labelbottom=True)

    # show x label for last man plot
    # axs_man[-1].tick_params(axis='x', labelbottom=True)
    # axs_man[-1].tick_params(axis='x', labelbottom=True)
    # ax.tick_params(axis='x', labelbottom=True)

    # set xlabel only for the last man plot
    # for ax in axs_man[:-1]:
    #    ax.xaxis.set_tick_params(labelbottom=False)
    #    # ax.tick_params(axis='x',labelbottom=False)
    #    ax.set_xlabel(None)

    if plot_gene:
        # plot nothing on axs[-2]
        axs[-2].axis('off')

        # TMP
        ax = axs[-1]
        gene_exon_ = format_gene_exon(gene_exon, chrompos, ax)
        print('xlim', ax.get_xlim())
        # plot gene
        plot_exon(ax, gene_exon_, figwidth)
        # plot_exon(ax, gene_exon_, figsize[0])
        # plot_exon(axs[-1],gene_exon_)
        # ax.tick_params(axis='both',labelbottom=False,labelleft=False)

        ax.axis('off')

    return fig


# for slide
# 1 row x several cols
def plot_locus_zoom_cols(ssd,
                         chrompos,
                         figtype,
                         phe_title,  # Phe
                         genetic_models=None,
                         legend=False,
                         plot_loss=False,
                         plot_gene=False,
                         annot=None,
                         gene_exon=None,
                         annot_pos=None,
                         ):

    # if len(chromposs) < 2:
    #    logger.error('chromposs {}'.format(chromposs))
    #    raise RuntimeError('chromposs<2:',len(chromposs))

    # print('boost plot nsnvs', len(ssobs_wgt[ssobs_wgt['boost_iteration'] != -1]))

    # print("boost_alpha", ssobs_wgt['boost_alpha'].head())

    # for figtype in [1,2]:
    # for figtype in [1,3,4]:
    # for figtype in ['1']:

    with sns.plotting_context('poster'):
        # +1 for loss
        gm_n = len(genetic_models) + 1
        # gm_n = 1

        fig, axs, figwidth = get_figsize_manhattan_plt_overlay_pos_multi(figtype, plot_gene, 1, gm_n, sharex=False, sharey=False)
        # fig, axs, figwidth = get_figsize_manhattan_plt_overlay_pos_multi(figtype, plot_gene, gm_n)

        # fig = plot_manhattan_multi(ssd, gene_exon, phe_title, chrompos, genetic_model, legend, plot_loss, plot_gene,
        fig = plot_manhattan_cols(ssd, gene_exon, phe_title, chrompos, genetic_models, legend, plot_loss, plot_gene,
                                  fig, axs, figwidth, annot, figtype, annot_pos=annot_pos)
        # TODO: add this?
        # ymax='0_default')

        # assume 2 axs for p-val manhattan
        man1_ylim = axs[1].get_ylim()
        man2_ylim = axs[2].get_ylim()
        man3_ylim = axs[3].get_ylim()
        # both ylim_min should be 0.0
        large_ylim_max = max(man1_ylim[1], man2_ylim[1], man3_ylim[1])
        # should be shared with other cols
        axs[1].set_ylim(ymax=large_ylim_max)
        axs[2].set_ylim(ymax=large_ylim_max)
        axs[3].set_ylim(ymax=large_ylim_max)

    # return (ylim_0, ylim_1)
    return fig


def plot_locus_zoom_multi(ssd,
                          chrompos,
                          figtype,
                          phe_title,  # Phe
                          genetic_models=None,
                          legend=False,
                          plot_loss=False,
                          plot_gene=False,
                          annot=None,
                          gene_exon=None,
                          ):

    # if len(chromposs) < 2:
    #    logger.error('chromposs {}'.format(chromposs))
    #    raise RuntimeError('chromposs<2:',len(chromposs))

    # print('boost plot nsnvs', len(ssobs_wgt[ssobs_wgt['boost_iteration'] != -1]))

    # print("boost_alpha", ssobs_wgt['boost_alpha'].head())

    # for figtype in [1,2]:
    # for figtype in [1,3,4]:
    # for figtype in ['1']:

    with sns.plotting_context('poster'):
        # +1 for loss
        gm_n = len(genetic_models) + 1

        fig, axs, figwidth = get_figsize_manhattan_plt_overlay_pos_multi(figtype, plot_gene, gm_n)

        fig = plot_manhattan_multi(ssd, gene_exon, phe_title, chrompos, genetic_models, legend, plot_loss, plot_gene,
                                   fig, axs, figwidth, annot, figtype)
        # TODO: add this?
        # ymax='0_default')

        # assume 2 axs for p-val manhattan
        man1_ylim = axs[1].get_ylim()
        man2_ylim = axs[2].get_ylim()
        man3_ylim = axs[3].get_ylim()
        # both ylim_min should be 0.0
        large_ylim_max = max(man1_ylim[1], man2_ylim[1], man3_ylim[1])
        # should be shared with other cols
        axs[1].set_ylim(ymax=large_ylim_max)
        axs[2].set_ylim(ymax=large_ylim_max)
        axs[3].set_ylim(ymax=large_ylim_max)

    # return (ylim_0, ylim_1)
    return fig


def plot_locus_zoom_multi_cols(ssd,
                               chromposs,
                               phe_title,  # Phe
                               legend=False,
                               plot_loss=False,
                               plot_gene=False,
                               gene_exon=None,
                               ylims=None):

    if len(chromposs) < 2:
        logger.error('chromposs {}'.format(chromposs))
        raise RuntimeError('chromposs<2:', len(chromposs), '. Use plot_locus_zoom_multi()')

    # print('boost plot nsnvs', len(ssobs_wgt[ssobs_wgt['boost_iteration'] != -1]))

    # print("boost_alpha", ssobs_wgt['boost_alpha'].head())

    # for figtype in [1,2]:
    # for figtype in [1,3,4]:
    for figtype in ['1']:

        with sns.plotting_context('poster'):

            # assume to have the same gm_n
            # +1 for loss
            gm_n = len(chromposs[0][1]) + 1
            ncols = len(chromposs)

            fig, axs, figwidth = get_figsize_manhattan_plt_overlay_pos_multi(figtype, plot_gene, gm_n, ncols)
            print(len(chromposs), len(axs[0]))
            assert (len(chromposs) == len(axs[0]))

            for i, chrompos_etc in enumerate(chromposs):
                chrompos = chrompos_etc[0]
                genetic_models = chrompos_etc[1]
                annot = chrompos_etc[2]
                axs_var = axs[:, i]
                fig = plot_manhattan_multi(ssd, gene_exon, phe_title, chrompos, genetic_models, legend, plot_loss, plot_gene,
                                           fig, axs_var, figwidth, annot, figtype, ymax='0_default')

            if ylims is None:
                # TODO: share y bet non-add p and add p
                # assume 2 axs for p-val manhattan

                # this doesn't seem to consider other ax even though ylim is shared
                large_ylim_max = 0.0
                for ax_coli in range(len(axs[0])):
                    man1_ylim = axs[1, ax_coli].get_ylim()
                    man2_ylim = axs[2, ax_coli].get_ylim()
                    large_ylim_max_coli = max(man1_ylim[1], man2_ylim[1])
                    large_ylim_max = max(large_ylim_max, large_ylim_max_coli)
                # man1_ylim=axs[1,0].get_ylim()
                # man2_ylim=axs[2,0].get_ylim()

                # for annotation height
                print('large_ylim_max', large_ylim_max)
                large_ylim_max = 1.1 * large_ylim_max
                print('large_ylim_max', large_ylim_max)

                # both ylim_min should be 0.0
                # large_ylim_max=max(man1_ylim[1],man2_ylim[1])
                # should be shared with other cols
                axs[1, 0].set_ylim(ymax=large_ylim_max)
                axs[2, 0].set_ylim(ymax=large_ylim_max)
            else:
                (ylim_0, ylim_1) = ylims
                axs[0, 0].set_ylim(ymax=ylim_0)
                axs[1, 0].set_ylim(ymax=ylim_1)
                axs[2, 0].set_ylim(ymax=ylim_1)

            # fname = fout + '.pval_alpha.size2.figtype' + str(figtype)
            # print('fname', fname)

            # fig.savefig(fname + '.bfr_layout.dpi1200.pdf', dpi=1200)

            # fig.subplots_adjust(hspace=0.05)
            # fig.subplots_adjust(hspace=0.07)
            # fig.subplots_adjust(hspace=0.15)
            # this makes hspace bet. man too large
            # fig.tight_layout()
            # fig.subplots_adjust(hspace=0.1)
            # fig.subplots_adjust(hspace=0.1,wspace=0.3)
            # fig.savefig(fname + '.png')
            # fig.savefig(fname + '.dpi1200.pdf', dpi=1200)
            # plt.close(fig)

            # ylim to return
            # ylim_0 = axs[0, 0].get_ylim()[1]
            # ylim_1 = axs[1, 0].get_ylim()[1]

    # return (ylim_0, ylim_1)
    return fig, axs
