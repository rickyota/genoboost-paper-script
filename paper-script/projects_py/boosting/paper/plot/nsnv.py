
import seaborn as sns
from logging import getLogger

# from ..validate.paper import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)



def plot_nsnv_box_strip_cat(
        stats,
        # plot_on,
        figtype=None,
        methods=None,  # to limit methods and order
        methodd=None,
        phes=None,  # to limit phes and order
        phed=None,
        # no_x=False,
        col_wrap=None,
        ylim=None,
        sharey=True,
        palette=None,
        palette_strip=None,
        markersize=None,
        ylog=False,
        ylabel=None,
        xtick_format=None,
        xtick_label_format=None,
        xtick_label=None,
        ytick=None,
        ytick_label=None,
        legend=None,
        legend_ncol=None,
        hlines=None,
        vlines=None,
):

    stats = stats.copy()

    if methods is None:
        methods = stats['method'].unique().tolist()

    stats.loc[:, 'Method'] = stats['method'].map(methodd)
    stats.loc[:, 'Phenotype'] = stats['phe'].map(phed)

    methodsp = pcommon.transit_plot(methods, methodd)
    phesp = pcommon.transit_plot(phes, phed)
    palettep = pcommon.transit_plot(palette, methodd)
    palettep_strip = pcommon.transit_plot(palette_strip, methodd)

    logger.debug('Method uniq {}'.format(stats['Method'].unique()))
    logger.debug('methodsp {}'.format(methodsp))

    kwargs_cat = dict(
        data=stats,
        x='Method',
        y='nsnv_use',
        # y='stat',
        # y=plot_on,
        order=methodsp,
        col='Phenotype',
        col_order=phesp,
        sharey=sharey,
        palette=palettep,
        col_wrap=col_wrap
    )

    # if no_x:
    #    pass
    #    # kwargs_cat.update(dict(x='Phenotype'))
    # else:
    # print('no_x!!')
    kwargs_cat |= dict(
        hue='Method',
        hue_order=methodsp,
        dodge=False
    )
    # kwargs_cat.update(dict(
    #    hue='Method',
    #    hue_order=methodsp,
    #    dodge=False
    # ))

    # if no_x:
    #    kwargs_cat.update(dict(x='Phenotype'))
    # else:
    #    kwargs_cat.update(dict(x='dummy', col='Phenotype', sharex=False,))

    _, fig_aspect = pcommon.parse_figtype(figtype)
    kwargs_cat |= fig_aspect

    # if col_wrap is not None:
    #    kwargs_cat |= dict(col_wrap=col_wrap)

    # logger.info('catplot')
    g = pcommon.plot_catplot_box_for_strip(**kwargs_cat)

    # somehow need this to show legend
    # should be here to show box legend
    g = pcommon.set_legend(g, legend, legend_ncol, methodsp)
    g = pcommon.change_box_color3(g)

    kwargs_swarm = dict(
        x='Method',
        y='nsnv_use',
        # y='stat',
        # y=plot_on,
        order=methodsp,
        # palette=['k', 'k'],
        # palette='dark:k',
        edgecolor='k',
        alpha=0.8,
    )

    if palettep_strip == 'k' or palettep_strip is None:
        kwargs_swarm |= dict(palette='dark:k')
        # TODO: make None default palette
        #elif palette_strip is not None:
    else:
        kwargs_swarm |= dict(palette=palettep_strip)

    # somehow warning
    # kwargs_swarm.update(dict(color='k'))
    # kwargs_swarm |= dict()
    # kwargs_swarm |= dict()

    # if no_x:
    #    pass
    # else:
    # no effect??
    kwargs_swarm |= dict(
        hue='Method',
        hue_order=methodsp,
        dodge=False
    )

    # needs if since default is size=5
    if markersize is not None:
        kwargs_swarm |= dict(size=markersize)

    # logger.info('swarm')
    g.map_dataframe(sns.swarmplot, **kwargs_swarm)

    g.set_titles(row_template='{row_name}', col_template='{col_name}')

    g = pcommon.set_ylog(g, ylog)

    g = pcommon.set_xtick(g, xtick_format)
    g = pcommon.set_xtick_label(g, xtick_label_format, xtick_label)
    g = pcommon.set_ytick(g, ytick, ytick_label)
    # g = pcommon.set_ytick_ntick(g, ytick_ntick)
    # g = set_ytick_suffix_auc(g, ytick_suffix)

    g = pcommon.set_xlabel(g, '')
    g = pcommon.set_ylabel(g, ylabel)

    # g = set_ylim_auc(g, ylim, plot_on)

    g = pcommon.set_hlines(g, hlines,
                           linestyle='-', linewidth=1.2, alpha=0.8,
                           color=pcommon.color_gray(),
                           zorder=0)

    g = pcommon.set_vlines(g, vlines,
                           linestyle='-', linewidth=1.2, alpha=0.5,
                           color=pcommon.color_gray(),
                           zorder=0)

    g = pcommon.set_move_legend(g, legend, legend_ncol, methods)

    # g.tight_layout()

    return g
