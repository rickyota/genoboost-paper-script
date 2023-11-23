

import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib as mpl
from logging import getLogger

# from ..validate.paper import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)


cb = sns.color_palette("colorblind", 10)
color_gray = cb[7]

# use **kwargs in this argument and pass directly to catplot()
# just indicate linewidth=1.0 and others
# TODO: This way is good, but how about indicating linewidth=1.0 in default argument?


def get_hue_order():
    # use [ori, add, lr] in this order here too
    sorter = [
        'GenoBoost',  # 'GenoBoost ori'
        'i. GenoBoost',  # 'GenoBoost ori'
        'A. GenoBoost',
        'Non-additive GenoBoost',
        'i. Non-additive GenoBoost',
        'Non-additive\nGenoBoost',
        'i. Non-additive\nGenoBoost',
        'ii. Non-additive\nGenoBoost',
        'a. Non-additive GenoBoost',
        'A. Non-additive GenoBoost',
        'a. Non-additive\nGenoBoost',
        'A. Non-additive\nGenoBoost',
        'Additive\neffects in\nnon-additive\nGenoBoost',
        'ii. Additive\neffects in\nnon-additive\nGenoBoost',
        'a1. Additive effects in\nnon-additive GenoBoost',
        'A1. Additive effects in\nnon-additive GenoBoost',
        'a1. Additive\neffects in\nnon-additive\nGenoBoost',
        'A1. Additive\neffects in\nnon-additive\nGenoBoost',
        'Dominance\neffects in\nnon-additive\nGenoBoost',
        'i. Dominance\neffects in\nnon-additive\nGenoBoost',
        'a2. Dominance effects in\nnon-additive GenoBoost',
        'A2. Dominance effects in\nnon-additive GenoBoost',
        'a2. Dominance\neffects in\nnon-additive\nGenoBoost',
        'A2. Dominance\neffects in\nnon-additive\nGenoBoost',
        'Additive GenoBoost',
        'ii. Additive GenoBoost',
        'Additive\nGenoBoost',
        'ii. Additive\nGenoBoost',
        'b. Additive GenoBoost',
        'B. Additive GenoBoost',
        'b. Additive\nGenoBoost',
        'B. Additive\nGenoBoost',
        'GenoBoost eps1',
        'GenoBoost eps2',
        'GenoBoost(add)',
        'GenoBoost(add) eps1',
        'GenoBoost(add) eps2',
        'GenoBoost(add) samepara',
        'GenoBoost(add) samepara eps1',
        'GenoBoost(lr)',
        'GenoBoost learning rate',
        'GenoBoost(lr) eps1',
        'LogitBoost (additive)',
        'snpnet',
        'ii. snpnet',
        'iii. snpnet',
        'snpnet (LASSO)',
        'ii. snpnet (LASSO)',
        'iii. snpnet (LASSO)',
        'snpnet\n(LASSO)',
        'ii. snpnet\n(LASSO)',
        'iii. snpnet\n(LASSO)',
        'B. snpnet\n(LASSO)',
        'LASSO',
        'iii. LASSO',
        'LASSO nocov',
        'iii. LASSO nocov',
        'LDAK-LASSO',
        'lassosum',
        'iii. lassosum',
        'iv. lassosum',
        'C. lassosum',
        'LDpred',
        'ldpred',
        'iv. LDpred',
        'v. LDpred',
        'D. LDpred',
        'PRS-CS',
        'prscs',
        'v. PRS-CS',
        'vi. PRS-CS',
        'PRScs',
        'vi. PRScs',
        'E. PRS-CS',
        'SBayesR',
        'vi. SBayesR',
        'vii. SBayesR',
        'vii. sbayesr',
        'F. SBayesR',
        'C+T',
        'vii. C+T',
        'viii. C+T',
        'G. C+T',
    ]
    return sorter


def get_hue_order_use(methods):
    # df['Methods']
    sorter = get_hue_order()
    # print('sorter', sorter)
    # print('methods', methods)
    for method in methods:
        if method not in sorter:
            raise RuntimeError('method not in sorter', method)
    sorter = [v for v in sorter if v in methods.tolist()]
    return sorter


def plot_catplot_box_for_strip_old(data, x, y, hue=None, col=None, row=None, kind=None,
                                   sharey=False, col_wrap=None, hue_order=None,
                                   height=None, aspect=None, legend=True, sharex=None, dodge=None, order=None, palette=None,
                                   facet_kws=None, linewidth=1.0, rc={}):
    # default whis
    kwargs = dict(x=x,
                  y=y,
                  hue=hue,
                  col=col,
                  row=row,
                  data=data,
                  kind='box',
                  sharey=sharey,
                  # linewidth=0,
                  linewidth=linewidth,  # thicker
                  # linewidth=1.0,  # thicker
                  # linewidth=2.5,  # thicker
                  fliersize=0.0,  # do not plot outlier
                  margin_titles=True,
                  col_wrap=col_wrap,
                  hue_order=hue_order,
                  legend=legend,
                  legend_out=False,  # False when use sns.move_legend()
                  # patch_artist=True, # sns already set this
                  )

    if height is not None:
        kwargs.update(dict(height=height))
    if aspect is not None:
        kwargs.update(dict(aspect=aspect))
    if sharex is not None:
        kwargs.update(dict(sharex=sharex))
    if dodge is not None:
        kwargs.update(dict(dodge=dodge))
    if order is not None:
        kwargs.update(dict(order=order))
    if palette is not None:
        kwargs.update(dict(palette=palette))
    if facet_kws is not None:
        kwargs.update(dict(facet_kws=facet_kws))

    print('rc', rc)
    print('kwargs', kwargs)
    g = sns.catplot(**rc, **kwargs)

    return g


def plot_catplot_box_for_strip(**kwargs):

    kwargs_default = dict(
        kind='box',
        sharey=False,
        linewidth=1.0,  # thicker
        fliersize=0.0,  # do not plot outlier
        margin_titles=True,
        legend_out=False,  # False when use sns.move_legend()
    )

    kwargs = kwargs_default | kwargs

    # print('kwargs', kwargs)
    g = sns.catplot(**kwargs)

    return g


def change_box_color3_ax(ax):
    # pathpatch = box
    # len(patchs) = x*y
    # len(ax.lines) = 6*x*y
    #
    # when have len(hue)=3, cat=2, color=c0,c1,c2
    # patchs: ??
    #

    # extract patchs of boxplot
    # https://stackoverflow.com/questions/72656861/how-to-add-hatches-to-boxplots-with-sns-boxplot-or-sns-catplot
    patchs = [patch for patch in ax.patches if isinstance(patch, mpl.patches.PathPatch)]
    # print('patches all', [type(patch) for patch in ax.patches])
    # print('patchs',patchs)

    for i, patch in enumerate(patchs):
        c = patch.get_facecolor()
        patch.set_edgecolor(c)

        # patch -> line
        k = 6 * i
        # for j in range(6*i,i*6+6):
        for l_ in [k, k + 1]:
            line = ax.lines[l_]
            # print("line", line)
            line.set_color(c)
            line.set_mfc(c)
            line.set_mec(c)
        for l_ in [k + 2, k + 3]:
            line = ax.lines[l_]
            # print("line", line)
            # transperent since setting white might overlap facecolor
            line.set_alpha(0)

    return ax


def change_box_color3(g):
    # make edgecolor and whisker vertical lines same as facecolor of box
    # only for g=catplot

    # print("g.axes", g.axes)
    for ax in g.axes.flatten():
        # print('ax title',ax.get_title())
        # print("ax", ax)
        # for sns version 0.11.2
        # https://matplotlib.org/stable/tutorials/intermediate/artists.html
        # print("ax artists", ax.artists)
        # print("ax artists", len(ax.artists))
        # print("ax lines", ax.lines)
        # print("ax patches", ax.patches)
        # print("ax patches len", len(ax.patches))

        ax = change_box_color3_ax(ax)

    return g


def set_xticklabels_newlines(ax):
    xticks = ax.get_xticklabels()
    print('xticks', xticks)

    def newline_xticks(xticks):
        # add newlines so that does not overwrap next ticks
        texts = [x.get_text() for x in xticks]
        print('texts', texts)

        def newlines_num(texts):
            def newline_in_text(text):
                return text.count('\n')

            newlines = []
            for i in range(len(texts)):
                if i % 2 == 0:
                    newlines += [0]
                elif i == len(texts) - 1:
                    num = newline_in_text(texts[i - 1])
                    newlines += [num + 1]
                else:
                    num_prev = newline_in_text(texts[i - 1])
                    num_next = newline_in_text(texts[i + 1])
                    newlines += [max(num_prev, num_next) + 1]
            return newlines

        newlines = newlines_num(texts)
        print('newlines', newlines)
        xticks = [''.join(['\n'] * n) + x for x, n in zip(texts, newlines)]
        return xticks
    xticks = newline_xticks(xticks)
    # xticks=[x.get_text() if i%2==0 else '\n\n'+x.get_text() for i,x in enumerate(xticks)]
    # xticks=[x.get_text() if i%2==0 else '\n\n'+x.get_text() for i,x in enumerate(xticks)]
    # xticks=[x.get_text() if i%2==0 else '\n\n'+x.get_text() for i,x in enumerate(xticks)]
    ax.set_xticklabels(xticks)
    return ax


def get_plot_on():
    return {'auc': 'AUC',
            'pauc_fpr0.1': 'pAUC',
            'pauc_fpr0.05': 'pAUC',
            'nagelkerke': 'Nagelkerke\'s $R^2$',
            'adj_mcfadden': 'Adjusted McFadden $R^2$'}


def get_legend_facetgrid(g):
    # check axes and find which is have legend
    leg = g.axes.flatten()[0].get_legend()
    if leg is not None:
        return leg
    # or legend may be on a figure
    # this happends when legend_out=True
    leg = g._legend
    if leg is None:
        raise RuntimeError
    return leg


def move_legend_facetgrid(g, **kwargs):
    sns.move_legend(g,
                    loc='upper center',
                    bbox_to_anchor=(0.5, 0.0),
                    title='Method',
                    # title=None,
                    # frameon=False,
                    frameon=True,
                    fancybox=False,
                    edgecolor='black',
                    handletextpad=0.3,
                    columnspacing=1.0,
                    # borderaxespad=0.0,
                    # mode='expand', # ng
                    **kwargs)
    # https://stackoverflow.com/questions/45201514/how-to-edit-a-seaborn-legend-title-and-labels-for-figure-level-functions

    print('facetgrid type', type(g))

    leg = get_legend_facetgrid(g)
    leg.get_frame().set_linewidth(1.0)
    return leg


def is_facetgrid(g):
    return isinstance(g, sns.FacetGrid)

# if type(ax)==sns.FacetGrid:


def get_legend_ax(ax):
    leg = None
    # check axes and find which has legend
    leg = ax.get_legend()
    if leg is not None:
        return leg
    if leg is None:
        raise NotImplementedError('uncomment code below')
    # not sure this could happen
    # or legend may be on a figure
    # leg = ax._legend
    return leg


def remove_legend_sns(g):
    if is_facetgrid(g):
        leg = get_legend_facetgrid(g)
    else:
        leg = get_legend_ax(g)
    leg.remove()


def plot_auc_cvs_phes_ulim_box_strip_cat(
        stats_auc,
        plot_on,
        modeld,
        figsize=(6, 6), xlim=None,
        ylim=None,
        logx=False,
        # xlabels_pos=None,
        col_wrap=None,
        methods=None,
        no_x=False, figtype=None,
        sharey=True,
        palette=None,
        ylabel=None,
        ysuffix=None,  # ytick_suffix
        linewidth=None,
        markersize=None,
        xtick=None,
        xtick_label=None,
        ytick_ntick=None,
        legend=None,
        legend_ncol=None,
        # legend_ncol=4
        legend_wrap=None,
        hlines=None,
        vline=None,
        phe_indent=True,
):

    stats_auc_ = stats_auc.copy()

    # TMP
    stats_auc_.loc[:, 'Method'] = stats_auc_['method']
    stats_auc_.loc[:, 'Phenotype'] = stats_auc_['phe']

    # stats_auc_.loc[:, 'Method'] = stats_auc_['method'].map(modeld)
    # if phe_indent:
    #    stats_auc_.loc[:, 'Phenotype'] = stats_auc_['phe'].map(get_phed_indent())
    # else:
    #    stats_auc_.loc[:, 'Phenotype'] = stats_auc_['phe'].map(get_phed())
    # stats_auc_.loc[:, 'Phenotype'] = stats_auc_['phe'].map(get_phed())
    # stats_auc_ = stats_auc_.sort_values(by='Method', key=method_sorter, kind='stable')
    stats_auc_.loc[:, 'dummy'] = 0

    if methods is None:
        methods_hue = None
    else:
        methods_hue = [modeld[x] for x in methods]

    print('Method uniq', stats_auc_['Method'].unique())
    print('order', get_hue_order_use(stats_auc_['Method']))

    # TMP: for legend test
    # g=sns.catplot(
    #    data=stats_auc_,
    #    x='Method',
    #    #order=get_hue_order_use(stats_auc_['Method']),
    #    y=plot_on,
    #    hue='Method',
    #    col='Phenotype',
    #    col_wrap=4,
    #    sharey=sharey,
    #    #sharey=True,
    #    palette=palette,
    #    kind='box',
    #    #legend=True,
    #    #legend_out=True,
    #    dodge=False,
    # )

    # somehow need this to show legend
    # g.add_legend()

    kwargs_cat = dict(
        data=stats_auc_,
        x='Method',
        order=get_hue_order_use(stats_auc_['Method']),
        y=plot_on,
        col='Phenotype',
        sharey=sharey,
        # sharey=True,
        palette=palette
    )

    if no_x:
        pass
        # kwargs_cat.update(dict(x='Phenotype'))
    else:
        print('no_x!!')
        kwargs_cat.update(dict(
            hue='Method',
            hue_order=get_hue_order_use(stats_auc_['Method']),
            dodge=False
        ))

    # if no_x:
    #    kwargs_cat.update(dict(x='Phenotype'))
    # else:
    #    kwargs_cat.update(dict(x='dummy', col='Phenotype', sharex=False,))

    if figtype == 1:
        # too large
        kwargs_cat.update(dict(height=8, aspect=1.5))
    elif figtype == 2:
        # for no_x
        kwargs_cat.update(dict(height=5, aspect=1.2))
        # kwargs_cat.update(dict(height=5, aspect=1))
    elif figtype == 3:
        # for not no_x
        kwargs_cat.update(dict(height=5, aspect=0.8))
    elif figtype == 4:
        kwargs_cat.update(dict(height=8, aspect=0.4))
    elif figtype == 5:
        kwargs_cat.update(dict(height=5, aspect=1.0))
    elif figtype == 6:
        kwargs_cat.update(dict(height=4, aspect=1.0))
        # kwargs_cat.update(dict(height=8, aspect=0.5))
    elif figtype == 7:
        kwargs_cat.update(dict(height=6, aspect=0.7))
    elif figtype == 8:
        kwargs_cat.update(dict(height=5, aspect=0.7))
    elif figtype == 9:
        kwargs_cat.update(dict(height=8, aspect=0.6))
    elif figtype == 10:
        kwargs_cat.update(dict(height=8, aspect=0.5))

    if linewidth is not None:
        kwargs_cat.update(dict(linewidth=linewidth))

    if col_wrap is not None:
        kwargs_cat.update(dict(col_wrap=col_wrap))

    g = plot_catplot_box_for_strip_old(**kwargs_cat,
                                       # somehow change facecolor too
                                       # rc=dict(boxprops=dict(color='w')), # to remove edge from legend
                                       )

    if legend is None:
        pass
    elif (legend == 'horizontal') and (legend_ncol is not None):
        # somehow need this to show legend
        # should be here to show box legend
        # https://stackoverflow.com/questions/28468584/seaborn-factorplot-set-series-order-of-display-in-legend
        # legend_ncol=4
        def get_legend_order_horizontal(labels, ncol):
            # get label[0], label[4], label[1], label[5]...
            labels_new = []
            nrow = (len(labels) + ncol - 1) // ncol
            # 0-4
            for i in range(ncol):
                # 0-1
                for j in range(nrow):
                    index = i + j * ncol
                    if index < len(labels):
                        labels_new.append(labels[index])
                    else:
                        break
            return labels_new
        legend_order = get_legend_order_horizontal(get_hue_order_use(stats_auc_['Method']), legend_ncol)
        g.add_legend(label_order=legend_order)
    else:
        g.add_legend()

    g = change_box_color3(g)

    kwargs = dict(
        x='Method',
        order=get_hue_order_use(stats_auc_['Method']),
        y=plot_on,
        alpha=0.8,
    )
    kwargs.update(dict(color='k'))

    if no_x:
        pass
    else:
        # no effect??
        kwargs.update(dict(
            hue='Method',
            hue_order=get_hue_order_use(stats_auc_['Method']),
            dodge=False
        ))

    # if no_x:
    #    kwargs.update(dict(x='Phenotype'))
    # else:
    #    kwargs.update(dict(x='dummy'))
    #    #kwargs.update(dict(x='dummy', col='Phenotype',))

    if markersize is not None:
        kwargs.update(dict(size=markersize))

    g.map_dataframe(sns.swarmplot, **kwargs)

    g.set_titles(row_template='{row_name}', col_template='{col_name}')

    # g.set_titles(col_template='{col_name}')
    # if not no_x:
    #    # set phenotype as x
    #    # g.set_titles(col_template='{col_name}')
    #    g.set_titles(col_template='')
    #    phes = stats_auc_['Phenotype'].unique().tolist()
    #    for ax, phe in zip(g.axes.flatten(), phes):
    #        print('get_xtick', ax.get_xticklabels())
    #        ax.set_xticklabels([phe])
    #        #ax.set(xticks=[0], xticklabels=[phe])

    if xtick is None:
        # access bottom axes only
        # for ax in grid._bottom_axes:
        for ax in g.axes.flatten():
            ax.xaxis.set_tick_params(rotation=90)
    elif xtick == 'norot':
        if xtick_label is not None:
            print('g.axes', g.axes.shape)
            for axis in g.axes.flatten():
                axis.tick_params(labelbottom=True)
            # if len(g.axes.shape) == 2:
            #    axes = g.axes[-1]
            # else:
            axes = g.axes
            # for ax in g.axes[-1]:
            for ax in axes:
                ax.set_xticklabels(xtick_label)
        else:
            pass
    elif xtick == 'newline':
        for ax in g.axes.flatten():
            set_xticklabels_newlines(ax)
    elif xtick == 'none':
        # should work
        # g.set_ticklabels(None)
        for ax in g.axes.flatten():
            ax.xaxis.set_tick_params(labelbottom=False)
            # ax.xaxis.set_tick_params(bottom=False,labelbottom=False)

    if ytick_ntick is not None:
        for ax in g.axes.flatten():
            ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick))

    if ylabel is None:
        ylabel = 'Stat'
        # ylabel = get_plot_on()[plot_on]

    g.set_xlabels('')
    # None does not remove label in g.set_xlabels()
    # g.set_xlabels(None)
    g.set_ylabels(ylabel)
    # g.set_axis_labels('', ylabel)
    # g.set_axis_labels('# SNVs', get_plot_on()[plot_on])
    # g.set_axis_labels('#SNVs', 'pAUC')

    if ysuffix is not None:
        def y_fmt(v, tick_number):
            # or round(v).is_integer()
            # if (100*v).is_integer():
            # should be all right for 3 digits
            # '5%' for integer and '5.1%' for float
            # print('v,100*v',v,100*v)
            if (round(100 * v, 3)).is_integer():
                # return '{:.0f}'.format(v) + '%
                return '{:.0%}'.format(v)
            else:
                return '{:.1%}'.format(v)
            # return str(v*100)+'%'
            # return '{:.1%}'.format(v)
            # return '{:2.1e}'.format(v)+'%'
            # return '{:2.1e}'.format(v)+'%'
            # return '{:2.2e}'.format(x).replace('e', 'x10^')

        for ax in g.axes.flatten():
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(y_fmt))

    # one snv only
    if ylim == 'bottom':
        # set minimum ylim
        if plot_on == 'accuracy':
            ylim_ = (0.90, 1)
        elif plot_on == 'auc' or plot_on.startswith('pauc'):
            ylim_ = (0.5, None)
        elif plot_on == 'nagelkerke':
            ylim_ = (0, 0, None)
        else:
            ylim_ = (0.0, None)
        for ax in g.axes.flatten():
            if ax.get_ylim()[0] < ylim_[0]:
                ax.set_ylim(ymin=ylim_[0])
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
                print('align0 ylim', ymin)
        else:
            # TODO
            print('Not Implemented')
    else:
        if ylim is None:
            if plot_on == 'accuracy':
                ylim_ = (0.90, 1)
            elif plot_on == 'auc' or plot_on.startswith('pauc'):
                ylim_ = (0.5, None)
            else:
                ylim_ = (0.0, None)
        elif ylim is False:
            ylim_ = None
        else:
            ylim_ = ylim

        if ylim_ is not False:
            # to avoid interrpting sharey
            g.set(ylim=ylim_)

    for ax in g.axes.flatten():
        if plot_on == 'nagelkerke':
            # TODO: set ticks
            pass
            # yticks = np.arange(0.0, 0.1, 0.005)
            # ymin, ymax = ax.get_ylim()
            # yticks = yticks[(yticks > ymin) & (yticks < ymax)]
            # print('ylim', ax.get_ylim())
            # print('yticks', yticks)
            # ax.set_yticks(yticks)  # minor=False
        else:
            # auc
            pass

    # print('sns.plotting_context()',sns.plotting_context())

    # make yticks fontsize smaller
    # for ax in g.axes.flatten():
    #    print('labelsize',ax.yaxis.get_tick_params(which='major'))
    #    fontsize=0.8*ax.yaxis.get_tick_params(which='major')['labelsize']
    #    ax.yaxis.set_tick_params(labelsize=fontsize)

    if hlines is not None:
        for hline_y in hlines:
            for ax in g.axes.flatten():
                ax.axhline(y=hline_y, linestyle='-', linewidth=1.2, alpha=0.8,
                           # color=sns.color_palette('colorblind')[0],
                           color=color_gray,
                           zorder=0)

    if vline is not None:
        for ax in g.axes.flatten():
            for vline_x in vline:
                ax.axvline(x=vline_x,
                           linestyle=':',
                           linewidth=1.2,
                           alpha=0.5,
                           color=color_gray,
                           zorder=0,  # move backward
                           )

    if legend is not None:
        kw = {}
        if legend == 'horizontal':
            if legend_ncol is not None:
                kw = dict(ncol=legend_ncol)
            else:
                # TODO: get from legend label
                # not working
                # if legend_out=True
                # legend_len = len(g.fig.legends[0].texts)
                # legend_len = len(g.fig.axes[0].legends.texts)
                legend_len = len(get_hue_order_use(stats_auc_['Method']))
                kw = dict(ncol=legend_len)
            kw.update(dict(borderaxespad=0.0))
            move_legend_facetgrid(g, **kw)
    else:
        remove_legend_sns(g)
        # sns.move_legend(g, loc='upper center', bbox_to_anchor=(0.5, 0.0),
        #                 title=None, frameon=False,
        #                 handletextpad=0.3,
        #                 columnspacing=1.0,
        #                 borderaxespad=0.0,
        #                #mode='expand', # ng
        #                 **kw)
        # sns.move_legend(g, loc='upper center', bbox_to_anchor=(0.5, 0.0),
        #                ncol=legend_ncol, title=None, frameon=False)

    # g.tight_layout()

    return g


def transit_plot(x, methodd):
    if methodd is None:
        return x
    else:
        if isinstance(x, dict):
            x = {methodd[k]: v for (k, v) in x.items()}
        elif isinstance(x, list):
            x = [methodd[k] for k in x]
        else:
            raise RuntimeError('Unknown type: ', type(x))
        return x


def parse_figtype(figtype):
    # (index)-(height)-(aspect)
    # figtype = '5-5-1.0'
    figindex = figtype.split('-')[0]
    height = float(figtype.split('-')[1])
    aspect = float(figtype.split('-')[2])
    return figindex, dict(height=height, aspect=aspect)


def set_xtick(g, xtick, xtick_label):
    if xtick is None:
        pass
    elif xtick == 'rot90':
        # access bottom axes only
        # for ax in grid._bottom_axes:
        for ax in g.axes.flatten():
            ax.xaxis.set_tick_params(rotation=90)
    elif xtick == 'norot':
        if xtick_label is None:
            pass
        else:
            print('g.axes', g.axes.shape)
            for axis in g.axes.flatten():
                axis.tick_params(labelbottom=True)
            # if len(g.axes.shape) == 2:
            #    axes = g.axes[-1]
            # else:
            axes = g.axes
            # for ax in g.axes[-1]:
            for ax in axes:
                ax.set_xticklabels(xtick_label)
    elif xtick == 'newline':
        for ax in g.axes.flatten():
            set_xticklabels_newlines(ax)
    elif xtick == 'none':
        # should work
        # g.set_ticklabels(None)
        for ax in g.axes.flatten():
            ax.xaxis.set_tick_params(labelbottom=False)
            # ax.xaxis.set_tick_params(bottom=False,labelbottom=False)
    return g


def set_ytick_ntick(g, ytick_ntick):
    if ytick_ntick is not None:
        for ax in g.axes.flatten():
            ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick))
    return g


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


def set_xlabel(g, xlabel):
    if xlabel is None:
        return g

    # if xlabel is False:
    #    xlabel = ''
    g.set_xlabels(xlabel)
    # None does not remove label in g.set_xlabels()
    # g.set_xlabels(None)
    return g


def set_ylabel(g, ylabel):
    if ylabel is None:
        return g

    # if ylabel is True:
    #    ylabel = 'Stat'
    g.set_ylabels(ylabel)
    return g


def set_margin_title_rot_font(g, margin_titile_font_prop):
    # rotate row titles
    # https://stackoverflow.com/questions/67427241/how-to-change-the-margin-title-color-in-a-seaborn-facetgrid
    # margin_title: plt.text.Annotation class https://matplotlib.org/3.1.1/api/text_api.html#matplotlib.text.Annotation
    # Annotation class inherits plt.text.Text class
    # Text class methods like .set_rotation() should be able to use
    for margin_title in g._margin_titles_texts:
        margin_title.set_rotation(0)
        # make the size twice
        fontsize = margin_title.get_fontsize()
        margin_title.set_fontsize(margin_titile_font_prop * fontsize)
        # if figtype in [1, 2, 3]:
        #    margin_title.set_fontsize(1.5 * fontsize)
        # else:
        #    margin_title.set_fontsize(2 * fontsize)
    return g


def set_bottom_label(g, blabel):
    if blabel is None:
        return g

    # if blabel is True:
    #    blabel = 'Stat'
    # g.set_ylabels(blabel)
    # g.set_ylabels(blabel)
    g.fig.suptitle(blabel, y=0.04)
    return g


def set_ylim_auc(g, ylim, plot_on):
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
def set_ytick_suffix(g, ytick_suffix):
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


def set_hlines(g, hlines, **kwargs):
    if hlines is not None:
        for hline_y in hlines:
            for ax in g.axes.flatten():
                ax.axhline(y=hline_y, **kwargs)
    return g


def set_vlines(g, vlines, **kwargs):
    if vlines is not None:
        for vline_x in vlines:
            for ax in g.axes.flatten():
                ax.axvline(x=vline_x, **kwargs)
    return g


def set_legend(g, legend, legend_ncol, methodsp):
    if legend is None:
        pass
    elif (legend == 'horizontal') and (legend_ncol is not None):
        # https://stackoverflow.com/questions/28468584/seaborn-factorplot-set-series-order-of-display-in-legend
        # legend_ncol=4
        def get_legend_order_horizontal(labels, ncol):
            # get label[0], label[4], label[1], label[5]...
            labels_new = []
            nrow = (len(labels) + ncol - 1) // ncol
            # 0-4
            for i in range(ncol):
                # 0-1
                for j in range(nrow):
                    index = i + j * ncol
                    if index < len(labels):
                        labels_new.append(labels[index])
                    else:
                        break
            return labels_new
        # legend_order = get_legend_order_horizontal(get_hue_order_use(stats['Method']), legend_ncol)
        legend_order = get_legend_order_horizontal(methodsp, legend_ncol)
        g.add_legend(label_order=legend_order)
    else:
        g.add_legend()
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
                           color=color_gray,
                           zorder=0)
            else:
                if plot_minor_hline:
                    ax.axhline(y=hline_y, linestyle='--', linewidth=0.8, alpha=0.5,
                               color=color_gray,
                               zorder=0)
    return g


def set_move_legend(g, legend, legend_ncol, methods):
    if legend is None:
        pass
    elif legend is False:
        remove_legend_sns(g)
    else:
        kw = {}
        if legend == 'horizontal':
            if legend_ncol is not None:
                kw = dict(ncol=legend_ncol)
            else:
                # TODO: get from legend label
                legend_len = len(methods)
                kw = dict(ncol=legend_len)
            kw.update(dict(borderaxespad=0.0))
            move_legend_facetgrid(g, **kw)
    return g


def plot_auc_box_strip_cat(
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

    methodsp = transit_plot(methods, methodd)
    phesp = transit_plot(phes, phed)
    palettep = transit_plot(palette, methodd)

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
        palette=palettep
    )

    # if no_x:
    #    pass
    #    # kwargs_cat.update(dict(x='Phenotype'))
    # else:
    # print('no_x!!')
    kwargs_cat |= dict(
        hue='Method',
        hue_order=methodsp,
        # should not necessary
        # dodge=False
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

    _, fig_aspect = parse_figtype(figtype)
    kwargs_cat |= fig_aspect

    if col_wrap is not None:
        kwargs_cat |= dict(col_wrap=col_wrap)

    # logger.info('catplot')
    g = plot_catplot_box_for_strip(**kwargs_cat)

    # somehow need this to show legend
    # should be here to show box legend
    g = set_legend(g, legend, legend_ncol, methodsp)
    g = change_box_color3(g)

    kwargs_swarm = dict(
        x='Method',
        y='stat',
        # y=plot_on,
        order=methodsp,
        palette=['k', 'k'],
        edgecolor='k',
        alpha=0.8,
    )
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

    if markersize is not None:
        kwargs_swarm |= dict(size=markersize)

    # logger.info('swarm')
    g.map_dataframe(sns.swarmplot, **kwargs_swarm)

    g.set_titles(row_template='{row_name}', col_template='{col_name}')

    g = set_xtick(g, xtick, xtick_label)
    g = set_ytick_ntick(g, ytick_ntick)
    g = set_ytick_suffix(g, ytick_suffix)

    g = set_xlabel(g, '')
    g = set_ylabel(g, ylabel)


    g = set_ylim_auc(g, ylim, plot_on)

    g = set_hlines(g, hlines,
                   linestyle='-', linewidth=1.2, alpha=0.8,
                   color=color_gray,
                   zorder=0)

    g = set_vlines(g, vlines,
                   linestyle='-', linewidth=1.2, alpha=0.5,
                   color=color_gray,
                   zorder=0)

    g = set_move_legend(g, legend, legend_ncol, methods)

    # g.tight_layout()

    return g


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
        margin_titile_font_prop=None
):

    stats = stats.copy()

    if methods is None:
        methods = stats['method'].unique().tolist()

    stats.loc[:, 'Method'] = stats['method'].map(methodd)
    stats.loc[:, 'Phenotype'] = stats['phe'].map(phed)
    stats.loc[:, 'Prop'] = stats['stat_on'].map(
        lambda x: x.replace('or_prop-', '') + '%'
    )
    propsp = [x + '%' for x in props]

    methodsp = transit_plot(methods, methodd)
    phesp = transit_plot(phes, phed)
    palettep = transit_plot(palette, methodd)

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

    _, fig_aspect = parse_figtype(figtype)
    kwargs_cat |= fig_aspect

    if col_wrap is not None:
        kwargs_cat |= dict(col_wrap=col_wrap)

    # logger.info('catplot')
    g = plot_catplot_box_for_strip(**kwargs_cat)

    # somehow need this to show legend
    # should be here to show box legend
    g = set_legend(g, legend, legend_ncol, methodsp)
    g = change_box_color3(g)

    kwargs_strip = dict(
        x='Prop',
        y='stat',
        order=propsp,
        hue='Method',
        hue_order=methodsp,
        palette=['k', 'k'],
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

    g = set_xlabel(g, xlabel)
    g = set_ylabel(g, ylabel)
    g = set_bottom_label(g, bottom_label)


    g = set_margin_title_rot_font(g, margin_titile_font_prop)

    g = set_hlines_odds_ratio(g, hlines)

    g = set_move_legend(g, legend, legend_ncol, methods)

    return g
