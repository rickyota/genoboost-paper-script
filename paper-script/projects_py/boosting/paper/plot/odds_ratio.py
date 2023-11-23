

import seaborn as sns
import matplotlib.ticker as ticker
from logging import getLogger

# from ..validate.paper import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)


def set_ticks_odds_ratio(g, ytick_ntick):
    # ytick_ntick = 4
    yticks = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
    for ax in g.axes.flatten():
        ymin, ymax = ax.get_ylim()
        yticks_ax = [v for v in yticks if ymin < v < ymax]
        # ok
        ax.yaxis.set_major_locator(ticker.FixedLocator(yticks_ax, nbins=ytick_ntick))
        # this makes sparser than ntick since counts tick outside of the lim
        # ax.yaxis.set_major_locator(ticker.FixedLocator(yticks,nbins=ytick_ntick))
        # ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick))
        # steps: allow 0.5 or 1.0 times
        # -> ng; this only allows
        # 1. [1.0, 2.0, 3.0,...]
        # 2. [0.5, 1.0, 1.5, 2.0, ...]
        # 3. x [1.0, 3.0, 5.0] -> since step=2
        # ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick,steps=[1,5,10]))
        # ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick,steps=[1,5,10]))
        # ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick))
    return g


def set_hlines_odds_ratio(g, hlines):
    for ax in g.axes.flatten():
        ymin, ymax = ax.get_ylim()
        if ymax - ymin < 4.5:
            plot_minor_hline = True
        else:
            plot_minor_hline = False
        for hline_y in hlines:
            if not (ymin < hline_y < ymax):
                continue
            # use list instead
            if hline_y in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]:
                ax.axhline(y=hline_y, linestyle='-', linewidth=1.2, alpha=0.8,
                           color=pcommon.color_gray(),
                           zorder=0)
            else:
                if plot_minor_hline:
                    ax.axhline(y=hline_y, linestyle='--', linewidth=0.8, alpha=0.5,
                               color=pcommon.color_gray(),
                               zorder=0)
    return g


def propp_format_int(x):
    x = x.replace('or_prop-', '')
    if float(x).is_integer():
        x = str(int(float(x)))
        # print('int',x)
        # error x='1.0'
        # x=str(int(x))
    # else:
        # print('float',x)
        # pass
    x = x + '%'
    return x


def plot_odds_ratio_box_strip_cat(
        stats,
        # plot_on,
        figtype=None,
        methods=None,  # to limit methods and order
        methodd=None,
        phes=None,  # to limit phes and order
        phed=None,
        props=None,
        # no_x=False,
        col_wrap=None,
        ylim=None,
        sharey=True,
        palette=None,
        markersize=None,
        strip_kind='swarm',
        xlabel=None,
        ylabel=None,
        bottom_label=None,
        xtick=None,
        xtick_label=None,
        ytick_ntick=None,
        ytick_suffix=None,
        legend=None,
        legend_ncol=None,
        hlines=None,
        vlines=None,
        # TODO: use plotting context instead
        margin_titile_font_prop=None
):

    stats = stats.copy()

    if methods is None:
        methods = stats['method'].unique().tolist()

    stats.loc[:, 'Method'] = stats['method'].map(methodd)
    stats.loc[:, 'Phenotype'] = stats['phe'].map(phed)

    stats.loc[:, 'Prop'] = stats['stat_on'].map(propp_format_int)
    propsp = [x + '%' for x in props]

    logger.debug('stats or: {}'.format(stats))
    # print('stats or', stats)

    methodsp = pcommon.transit_plot(methods, methodd)
    phesp = pcommon.transit_plot(phes, phed)
    palettep = pcommon.transit_plot(palette, methodd)

    kwargs_cat = dict(
        data=stats,
        x='Prop',
        y='stat',
        # y=plot_on,
        order=propsp,
        hue='Method',
        hue_order=methodsp,
        col='Phenotype',
        col_order=phesp,
        sharey=sharey,
        palette=palettep,
        # dodge=False
    )

    figindex, fig_aspect = pcommon.parse_figtype(figtype)
    kwargs_cat |= fig_aspect

    if col_wrap is not None:
        kwargs_cat |= dict(col_wrap=col_wrap)

    # logger.info('catplot')
    g = pcommon.plot_catplot_box_for_strip(**kwargs_cat)

    # somehow need this to show legend
    # should be here to show box legend
    g = pcommon.set_legend(g, legend, legend_ncol, methodsp)
    g = pcommon.change_box_color3(g)

    kwargs_strip = dict(
        x='Prop',
        y='stat',
        order=propsp,
        hue='Method',
        hue_order=methodsp,
        palette='dark:k',
        edgecolor='k',
        alpha=0.8,
        dodge=True
    )

    if markersize is not None:
        kwargs_strip |= dict(size=markersize)

    if strip_kind == 'swarm':
        g.map_dataframe(sns.swarmplot, **kwargs_strip)
    elif strip_kind == 'strip':
        kwargs_strip |= dict(jitter=0)
        g.map_dataframe(sns.stripplot, **kwargs_strip)
    else:
        raise RuntimeError('Unknown strip_kind', strip_kind)

    g.set_titles(row_template='{row_name}', col_template='{col_name}')

    g = set_ticks_odds_ratio(g, ytick_ntick)

    g = pcommon.set_xlabel(g, xlabel)
    g = pcommon.set_ylabel(g, ylabel)
    if figindex == '9':
        g = pcommon.set_bottom_label(g, bottom_label, 0.02)
    else:
        g = pcommon.set_bottom_label(g, bottom_label, 0.04)

    g = pcommon.set_margin_title_rot_font(g, margin_titile_font_prop)

    g = set_hlines_odds_ratio(g, hlines)

    g = pcommon.set_move_legend(g, legend, legend_ncol, methods)

    return g
