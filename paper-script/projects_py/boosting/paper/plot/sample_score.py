
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from statannotations.Annotator import Annotator
import numpy as np
from logging import getLogger

# from ..validate.paper import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)


def get_figsize_sample_summary_both(figtype, summary, mode=None):
    nrows = len(summary)
    kws = {}
    if nrows == 3:
        # assume to have box
        if figtype == '1':
            # figsize = (5,10)
            figsize = (8, 16)
            height_ratios = [5, 2, 5]
        elif figtype == '2':
            figsize = (12, 12)
            # expand x for legend
            # figsize = (8.5,12)
            # figsize = (8,12)
            # figsize = (9,12)
            height_ratios = [5, 2, 5]
        elif figtype == '3':
            figsize = (12, 10)
            # expand x for legend
            height_ratios = [5, 2, 5]
        else:
            raise NotImplementedError('figtype: ', figtype)
        kws |= dict(gridspec_kw={'height_ratios': height_ratios})
        # kws|=('gridspec_kw',{'height_ratios': height_ratios})

    fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=figsize,
                            sharex='col',
                            sharey='row',
                            # layout="constrained",  # ng
                            squeeze=False,
                            **kws
                            )

    return fig, axs


def axplot_sample_hist(ax, rs, stat, x, xlabel, palette=None, mode=None, legend=True):
    if stat == 'proportion':
        ylabel = 'Frequency of samples'
    elif stat == 'count':
        ylabel = 'Count of samples'
    else:
        ylabel = stat

    # palette = paletted_diesase_status()

    #logger.debug('rs {}'.format(rs))
    #print('inf in ?',(~np.isfinite(rs[x])).any())
    #print('inf in ?',rs.loc[~np.isfinite(rs[x]),x])

    ax = sns.histplot(data=rs,
                      x=x,
                      # x='score',
                      hue='Phenotype',
                      hue_order=['Case', 'Control'],
                      palette=palette,
                      stat=stat,
                      common_norm=False,
                      ax=ax,
                      edgecolor='w',
                      #linewidth=0.0
                      )

    if mode == 'simple':
        ax.set_xlabel('GenoBoost sample score')
    else:
        ax.set_xlabel(xlabel)

    # xlim = (rs['score'].min() * 1.05, rs['score'].max() * 1.05)
    # ax.set_xlim(xlim)
    ax.set_ylabel(ylabel)

    if mode == 'simple':
        ax.yaxis.set_tick_params(which='both', left=False, labelleft=False)

    if mode == 'simple' or legend is False:
        pcommon.set_move_legend_ax(ax, False)
        # remove_legend_sns(ax)
    else:
        pcommon.set_move_legend_ax(ax, True)
        # move_legend_ax(ax)
    # sns.move_legend(ax)

    if mode == 'simple':
        texty = ax.get_ylim()[1] * 0.7
        ax.text(-1.0, texty, 'Control', color=palette['Control'], horizontalalignment='center',)
        ax.text(1.0, texty, 'Case', color=palette['Case'], horizontalalignment='center',)

    return ax


def axplot_sample_box(ax, rs, x, xlabel, palette):
    # x = 'score'
    y = 'Phenotype'
    order = ['Case', 'Control']
    # order = ['Control', 'Case']

    # palette = paletted_diesase_status()
    ax = sns.boxplot(x=x, y=y, order=order, data=rs,
                     palette=palette, ax=ax)

    pairs = [('Control', 'Case')]
    annotator = Annotator(ax, pairs, data=rs, x=x, y=y, order=order, orient='h')
    # annotator = Annotator(ax, pairs, data=rs, x=x, y=y)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
    annotator.apply_and_annotate()

    ax.set_xlabel(xlabel)
    ax.set_ylabel(None)

    return ax


def axplot_sample_kde(ax, rs, x, y, palette, xlabel=None, ylabel=None, legend_kw={},):
    print('x,y', x, y)
    # palette = paletted_diesase_status()

    # takes time
    ax = sns.kdeplot(
        data=rs,
        x=x,
        y=y,
        hue='Phenotype',
        hue_order=['Case', 'Control'],
        common_norm=False,
        # TODO
        levels=[0.1] + list(np.linspace(0.2, 1.0, 5)),
        # levels=[0.02, 0.05, 0.1] + list(np.linspace(0.2, 1.0, 5)),
        # levels=[0.05] + list(np.linspace(0.1, 1.0, 5)),
        # levels=[0.02, 0.06, ] + list(np.linspace(0.1, 1.0, 5)),
        # levels=[0.02, 0.04, 0.06, 0.08] + list(np.linspace(0.1, 1.0, 10)),
        # thresh=0.001,
        # thresh=0.0,
        palette=palette,
        linewidths=2.0,  # change here so that not chaing legend linewidth
        ax=ax,
    )

    cb = sns.color_palette("colorblind", 10)
    ax = sns.regplot(
        data=rs,
        x=x,
        y=y,
        ci=None,
        scatter=False,
        color=cb[2],
        line_kws={'linewidth': 1.8},
        ax=ax
    )

    # ax.get_legend().set_title(None)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    pcommon.move_legend_ax(ax,
                           # loc="upper left",
                           # bbox_to_anchor=(1.05, 1),
                           **legend_kw
                           )

    return ax


def axplot_prop_vline(ax, rs, x, prop):
    # vertical line at top n proportion
    n_topprop = int(len(rs) * prop)
    print('n_topprop', n_topprop)
    score_topprop = rs[x].sort_values(ascending=False).iloc[n_topprop]
    print('score_topprop', score_topprop)

    # for prop
    # ax.plot([score_topprop, score_topprop], [0.0, 0.01], linewidth=0.3)
    ax.axvline(x=score_topprop)
    return ax


# def validate_sample_summary_both(fout, rs, x, y,
def plot_sample_score_summary(rs, x, y,
                              phenod=None,
                              xlabel=None, ylabel=None,
                              summary=['hist', 'box', 'kde'],
                              palette=None,
                              mode=None,
                              # phenotype=None,
                              title=None,
                              plot_prop=False,
                              stat=None,
                              figtype=None,
                              ):

    sharex_rs = True

    palettep = pcommon.transit_plot(palette, phenod)

    # if phenotype is None:
    #    title = ''
    # else:
    #    title = str(get_phed()[phenotype])

    # for stat in ['count', 'proportion']:

    # for figtype in [2]:
    # for figtype in [3]:

    # with sns.plotting_context('poster',
    #                          #rc={'lines.linewidth': 2.0}
    #                          # rc={'ytick.labelsize': 18, 'legend.fontsize':24}
    #                          ):

    rs.loc[:, 'Phenotype'] = rs['status'].map(phenod)

    fig, axs = get_figsize_sample_summary_both(figtype, summary, mode)

    axs_0 = axs[:, 0]
    for ax, plotkind in zip(axs_0, summary):
        if plotkind == 'hist':
            axplot_sample_hist(ax, rs, stat, x, xlabel, palettep, mode, legend=False)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
        elif plotkind == 'box':
            axplot_sample_box(ax, rs, x, xlabel, palettep)
        elif plotkind == 'kde':
            axplot_sample_kde(ax, rs, x, y, palettep, xlabel, ylabel, legend_kw=dict(loc='lower left', bbox_to_anchor=(1.05, 0)))
            # make these mode?
            ax.xaxis.set_major_locator(ticker.MaxNLocator(3))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
            ax.xaxis.set_tick_params(which='minor', bottom=False)
            ax.yaxis.set_tick_params(which='minor', bottom=False)

    if plot_prop is not False:
        prop = plot_prop
        for ax in axs_0:
            axplot_prop_vline(ax, rs, x, prop)

    # plot y
    assert summary[2] == 'kde'
    axs_1 = axs[:-1, 1]
    # axs_1=axs[:,1]
    for ax, plotkind in zip(axs_1, summary):
        if plotkind == 'hist':
            axplot_sample_hist(ax, rs, stat, y, ylabel, palettep, mode, legend=True)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
        elif plotkind == 'box':
            axplot_sample_box(ax, rs, y, ylabel, palettep)

    if plot_prop is not False:
        prop = plot_prop
        for ax in axs_1:
            axplot_prop_vline(ax, rs, y, prop)

    axs[-1, 1].axis('off')

    if sharex_rs:
        xlim0 = axs[0, 0].get_xlim()
        xlim1 = axs[0, 1].get_xlim()
        xlim = (min(xlim0[0], xlim1[0]), max(xlim0[1], xlim1[1]))
        for ax in axs.flatten():
            ax.set_xlim(xlim)
        # set ylim for kde
        axs[2, 1].set_ylim(xlim)

    # remove x label
    # if nrows>=2:
    for ax in axs[0, :-1]:
        ax.set_xlabel(None)
    axs[1, 0].set_xlabel(None)
    for ax in axs[1, :]:
        ax.set_ylabel(None)

    if title is not None:
        # axs[0].set_title(title)
        fig.suptitle(title)

    # equal axis and square box
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axes_box_aspect.html
    # this makes axis skewed
    # axs[-1].axis('equal')
    # this only changes axs[-1]
    # axs[-1].set_box_aspect(1)
    # axs[-1].set_aspect('equal', 'box')

    return fig
