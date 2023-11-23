

import seaborn as sns
import matplotlib.ticker as ticker
from logging import getLogger

# from ..validate.paper import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)


# moved to config
# def get_plot_on():
#    return {'auc': 'AUC',
#            'pauc_fpr0.1': 'pAUC',
#            'pauc_fpr0.05': 'pAUC',
#            'nagelkerke': 'Nagelkerke\'s $R^2$',
#            'adj_mcfadden': 'Adjusted McFadden $R^2$'}


def set_ylim_acc(g, ylim, plot_on):
    if ylim is None:
        pass
    elif ylim == 'bottom' or ylim == 'fix':
        # set minimum ylim
        if plot_on == 'accuracy':
            ylim_ = (0.90, 1)
        elif plot_on == 'auc' or plot_on.startswith('pauc'):
            ylim_ = (0.5, None)
        elif plot_on == 'nagelkerke':
            ylim_ = (0.0, None)
        else:
            ylim_ = (0.0, None)

        if ylim == 'fix':
            g.set(ylim=ylim_)
        elif ylim == 'bottom':
            for ax in g.axes.flatten():
                if ax.get_ylim()[0] < ylim_[0]:
                    ax.set_ylim(ymin=ylim_[0])
        else:
            raise NotImplementedError
    elif ylim == 'align0':
        if plot_on == 'nagelkerke':
            # make length of <0.0 of y axes the same
            # this can align y=0.0
            ymin_max = -0.002
            ymax_max = max([ax.get_ylim()[1] for ax in g.axes.flatten()])
            # for simplicity
            prop_neg = ymin_max / ymax_max
            # prop_neg=-ymin_max/(ymax_max-ymin_max)
            for ax in g.axes.flatten():
                ymax = ax.get_ylim()[1]
                ymin = ymax * prop_neg
                ax.set_ylim(ymin=ymin)
                # print('align0 ylim', ymin)
        else:
            # TODO
            logger.info('Not Implemented')
    else:
        g.set(ylim=ylim)
    return g


# TODO: if ytick_suffix=='%
def set_ytick_suffix_acc(g, ytick_suffix):
    if ytick_suffix is not None:
        def y_fmt(v, tick_number):
            # or round(v).is_integer()
            # if (100*v).is_integer():
            # should be all right for 3 digits
            # '5%' for integer and '5.1%' for float
            # print('v,100*v',v,100*v)
            if (round(100 * v, 3)).is_integer():
                return '{:.0%}'.format(v)
            else:
                return '{:.1%}'.format(v)

        for ax in g.axes.flatten():
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_fmt))
    return g


# exclude None if it does not work without the argument
def plot_acc_box_strip_cat(
        stats,
        # plot_on,
        plot_on,
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
        markersize=None,
        ylabel=None,
        xtick_format=None,
        xtick_label_format=None,
        xtick_label=None,
        ytick_ntick=None,
        ytick_suffix=None,
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
    # stats.loc[:, 'Phenotype'] = stats['phe'].map(phed)

    methodsp = pcommon.transit_plot(methods, methodd)
    phesp = pcommon.transit_plot(phes, phed)
    palettep = pcommon.transit_plot(palette, methodd)

    logger.debug('Method uniq {}'.format(stats['Method'].unique()))
    logger.debug('methodsp {}'.format(methodsp))

    kwargs_cat = dict(
        data=stats,
        x='Method',
        y='stat',
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
        y='stat',
        # y=plot_on,
        order=methodsp,
        palette='dark:k',
        edgecolor='k',
        alpha=0.8,
    )
    # somehow warning -> palette='dark:k'
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

    g = pcommon.set_xtick(g, xtick_format)
    g = pcommon.set_xtick_label(g, xtick_label_format, xtick_label)
    g = pcommon.set_ytick_ntick(g, ytick_ntick)
    g = set_ytick_suffix_acc(g, ytick_suffix)

    g = pcommon.set_xlabel(g, '')
    g = pcommon.set_ylabel(g, ylabel)

    g = set_ylim_acc(g, ylim, plot_on)

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


def plot_acc_simu_box_strip_cat(
        stats,
        # plot_on,
        plot_on,
        figtype=None,
        methods=None,  # to limit methods and order
        methodd=None,
        phes=None,  # to limit phes and order
        phed=None,
        # no_x=False,
        # col_wrap=None,
        ylim=None,
        sharey=True,
        palette=None,
        markersize=None,
        ylabel=None,
        xtick=None,
        xtick_label=None,
        ytick_ntick=None,
        ytick_suffix=None,
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
    # stats.loc[:, 'Phenotype'] = stats['phe'].map(phed)

    methodsp = pcommon.transit_plot(methods, methodd)
    phesp = pcommon.transit_plot(phes, phed)
    palettep = pcommon.transit_plot(palette, methodd)

    logger.debug('Method uniq {}'.format(stats['Method'].unique()))
    logger.debug('methodsp {}'.format(methodsp))

    kwargs_cat = dict(
        data=stats,
        x='Method',
        y='stat',
        # y=plot_on,
        order=methodsp,
        col='h2dom',
        # col='Phenotype',
        # col_order=phesp,
        # col_order=['0.05', '0'],
        row='h2',
        sharey=sharey,
        # palette=palettep,
        # col_wrap=col_wrap
    )

    # if no_x:
    #    pass
    #    # kwargs_cat.update(dict(x='Phenotype'))
    # else:
    # print('no_x!!')
    kwargs_cat |= dict(
        hue='ncausal',
        hue_order=['100', '1000'],
        # hue='Method',
        # hue_order=methodsp,
        dodge=True
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
    # g.add_legend()
    # g = pcommon.set_legend(g, legend, legend_ncol, ['100', '1000'])
    g = pcommon.set_legend(g, legend, legend_ncol, methodsp)
    g = pcommon.change_box_color3(g)

    # g.add_legend(title='# causal variants')

    kwargs_swarm = dict(
        x='Method',
        y='stat',
        # y=plot_on,
        order=methodsp,
        palette='dark:k',
        edgecolor='k',
        alpha=0.8,
    )
    # somehow warning -> palette='dark:k'
    # kwargs_swarm.update(dict(color='k'))
    # kwargs_swarm |= dict()
    # kwargs_swarm |= dict()

    # if no_x:
    #    pass
    # else:
    # no effect??
    kwargs_swarm |= dict(
        hue='ncausal',
        hue_order=['100', '1000'],
        # hue='Method',
        # hue_order=methodsp,
        dodge=True,
        # dodge=False
        # no effect
        # legend=True,
    )

    # needs if since default is size=5
    if markersize is not None:
        kwargs_swarm |= dict(size=markersize)

    # logger.info('swarm')
    g.map_dataframe(sns.swarmplot, **kwargs_swarm)

    g.set_titles(row_template='$h^2$={row_name}', col_template='$h^2_d$={col_name}')
    # g.set_titles(row_template='Additive + Dominance $h^2$={row_name}', col_template='Dominance $h^2$={col_name}')
    # g.set_titles(row_template='$h^2$={row_name}', col_template='$Dominance h^2_{dom}$={col_name}')

    g = pcommon.set_xtick_label(g, xtick, xtick_label)
    g = pcommon.set_ytick_ntick(g, ytick_ntick)
    g = set_ytick_suffix_acc(g, ytick_suffix)

    g = pcommon.set_xlabel(g, '')
    g = pcommon.set_ylabel(g, ylabel)

    g = set_ylim_acc(g, ylim, plot_on)

    g = pcommon.set_hlines(g, hlines,
                           linestyle='--', linewidth=1.2, alpha=0.8,
                           color=pcommon.color_gray(),
                           zorder=0)

    g = pcommon.set_vlines(g, vlines,
                           linestyle='-', linewidth=1.2, alpha=0.5,
                           color=pcommon.color_gray(),
                           zorder=0)
    # g.add_legend()

    # TMP commnet
    g = pcommon.set_move_legend(g, legend, legend_ncol, ['100', '1000'], title='# causal variants',
                                place=None,
                                loc='lower left',
                                bbox_to_anchor=(1.05, 0.15))
    # place='lower right outside')
    # g = pcommon.set_move_legend(g, legend, legend_ncol, methods)

    # g.tight_layout()

    return g
