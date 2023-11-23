

import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib as mpl
from logging import getLogger

# from ..validate.paper import config as pconfig

logger = getLogger(__name__)


def color_gray():
    cb = sns.color_palette("colorblind", 10)
    color_gray = cb[7]
    return color_gray


def transit_plot(x, methodd):
    if methodd is None:
        return x
    else:
        if x is None:
            return x
        elif isinstance(x, str):
            return x
        elif isinstance(x, dict):
            # FIXME: raise eror message when k is not in methodd
            x = {methodd[k]: v for (k, v) in x.items() if k in methodd}
            # x = {methodd[k]: v for (k, v) in x.items()}
        elif isinstance(x, list):
            # not tested
            x = [methodd[k] for k in x if k in methodd]
            # x = [methodd[k] for k in x]
        else:
            raise RuntimeError('Unknown type: ', type(x))
        return x


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


def _move_legend_facetgrid(g, place, **kwargs):
    if place is None:
        pass
    elif place == 'lower center outside':
        kwargs |= dict(loc='upper center',
                       bbox_to_anchor=(0.5, 0.0))
    elif place == 'lower right outside':
        kwargs |= dict(loc='lower left',
                       bbox_to_anchor=(1.05, 0.0))
    else:
        raise RuntimeError('Unknown place: {}'.format(place))

    sns.move_legend(g,
                    # loc='upper center',
                    # bbox_to_anchor=(0.5, 0.0),
                    # title should be in kwargs
                    # title='Method',
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

    # logger.debug('facetgrid type {}'.format(type(g)))

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


# def transit_plot(x, methodd):
#    if methodd is None:
#        return x
#    else:
#        if isinstance(x, dict):
#            x = {methodd[k]: v for (k, v) in x.items()}
#        elif isinstance(x, list):
#            x = [methodd[k] for k in x]
#        else:
#            raise RuntimeError('Unknown type: ', type(x))
#        return x


def parse_figtype(figtype, return_figsize=False):
    # (index)-(height)-(aspect)
    # figtype = '5-5-1.0'
    figindex = figtype.split('-')[0]
    height = float(figtype.split('-')[1])
    aspect = float(figtype.split('-')[2])
    if return_figsize:
        # for matplotlib
        return figindex, (height * aspect, height)
    else:
        # for seaborn
        return figindex, dict(height=height, aspect=aspect)


def set_xticklabels_newlines(ax):
    xticks = ax.get_xticklabels()
    logger.debug('xticks {}'.format(xticks))

    def newline_xticks(xticks):
        # add newlines so that does not overwrap next ticks
        texts = [x.get_text() for x in xticks]
        logger.debug('texts {}'.format(texts))

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
        logger.debug('newlines {}'.format(newlines))
        xticks = [''.join(['\n'] * n) + x for x, n in zip(texts, newlines)]
        return xticks
    xticks = newline_xticks(xticks)
    # xticks=[x.get_text() if i%2==0 else '\n\n'+x.get_text() for i,x in enumerate(xticks)]
    # xticks=[x.get_text() if i%2==0 else '\n\n'+x.get_text() for i,x in enumerate(xticks)]
    # xticks=[x.get_text() if i%2==0 else '\n\n'+x.get_text() for i,x in enumerate(xticks)]
    ax.set_xticklabels(xticks)
    return ax


# xtick and its label
# format -> type? method?
def set_xtick_label(g, xtick_label_format, xtick_label):
    if xtick_label_format is None:
        pass
    elif xtick_label_format == 'rot90':
        # access bottom axes only
        # for ax in grid._bottom_axes:
        for ax in g.axes.flatten():
            ax.xaxis.set_tick_params(rotation=90)
    elif xtick_label_format == 'norot':
        if xtick_label is None:
            pass
        else:
            # logger.debug('g.axes {}'.format(g.axes.shape))
            for axis in g.axes.flatten():
                axis.tick_params(labelbottom=True)
            # if len(g.axes.shape) == 2:
            #    axes = g.axes[-1]
            # else:
            # axes = g.axes
            # TODO: ok?
            axes = g.axes.flatten()
            # for ax in g.axes[-1]:
            for ax in axes:
                ax.set_xticklabels(xtick_label)
    elif xtick_label_format == 'newline':
        for ax in g.axes.flatten():
            set_xticklabels_newlines(ax)
    elif xtick_label_format == 'none':
        # should work
        # g.set_ticklabels(None)
        for ax in g.axes.flatten():
            ax.xaxis.set_tick_params(labelbottom=False)
            # ax.xaxis.set_tick_params(bottom=False,labelbottom=False)
    return g


def set_xtick(g, xtick_format):
    if xtick_format is None:
        pass
    elif xtick_format == 'none':
        # should work
        # g.set_ticklabels(None)
        for ax in g.axes.flatten():
            # this  also remove label?
            ax.xaxis.set_tick_params(bottom=False)
            # ax.xaxis.set_tick_params(bottom=False,labelbottom=False)
    return g


def set_ytick_ntick(g, ytick_ntick):
    if ytick_ntick is not None:
        for ax in g.axes.flatten():
            ax.yaxis.set_major_locator(ticker.MaxNLocator(ytick_ntick))
    return g


# TODO: rename to set_ytick_label()
def set_ytick(g, ytick, ytick_label):
    if ytick is not None:
        for ax in g.axes.flatten():
            ax.set_yticks(ytick)

    if ytick_label is not None:
        for ax in g.axes.flatten():
            ax.set_yticklabels(ytick_label)
    # for ax in g.axes.flatten():
    #    ax.xaxis.set_tick_params(rotation=90)
    # remove x ticks
    # for ax in g.axes.flatten():
    #    ax.axes.xaxis.set_ticks([])

    # TODO: arg
    # remove y minor ticks
    for ax in g.axes.flatten():
        # ax.minorticks_off()
        ax.yaxis.set_minor_locator(ticker.NullLocator())

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


def set_ylog(g, ylog):
    if ylog:
        # for ax in g.axes.flatten():
        # ax = g.axes.flatten()[0]
        for ax in g.axes.flatten():
            ax.set_yscale('log')
        # error
        # g.set_yscale('log')
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


def set_bottom_label(g, label, y):
    if label is None:
        return g

    # if blabel is True:
    #    blabel = 'Stat'
    # g.set_ylabels(blabel)
    # g.set_ylabels(blabel)
    # g.fig.suptitle(blabel, y=0.04)
    g.fig.suptitle(label, y=y)
    return g


def set_hlines(g, hlines, **kwargs):
    if hlines is not None:
        for ax in g.axes.flatten():
            ymin, ymax = ax.get_ylim()
            for hline_y in hlines:
                if not (ymin < hline_y < ymax):
                    continue
                ax.axhline(y=hline_y, **kwargs)
    return g


def set_vlines(g, vlines, **kwargs):
    if vlines is not None:
        for ax in g.axes.flatten():
            for vline_x in vlines:
                ax.axvline(x=vline_x, **kwargs)
    return g


def set_legend(g, legend, legend_ncol, methodsp):
    if legend is None:
        # None means doing nth
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
        # no effect
        # g.add_legend(label_order=legend_order, titile='# causal variants')
    else:
        # if legend=False, will be removed in set_move_legend()
        g.add_legend()
    return g


# TODO: remove title='Method' and rewrite other fn
# FIXME: change default place=None, otherwise legend will be renewed
def set_move_legend(g, legend, legend_ncol, methods, title='Method', place='lower center outside', **kwargs):
    if legend is None:
        pass
    elif legend is False:
        remove_legend_sns(g)
    else:
        # kw = {}
        kw = kwargs
        if legend == 'horizontal':
            if legend_ncol is not None:
                kw |= dict(ncol=legend_ncol)
            else:
                # TODO: get from legend label
                legend_len = len(methods)
                kw |= dict(ncol=legend_len)
            kw.update(dict(borderaxespad=0.0))
            if title is not None:
                kw.update(dict(title=title))
            _move_legend_facetgrid(g, place, **kw)
    return g

# def set_move_legend(g, legend, legend_ncol, methods, title='Method'):
#    if legend is None:
#        pass
#    elif legend is False:
#        remove_legend_sns(g)
#    else:
#        kw = {}
#        if legend == 'horizontal':
#            if legend_ncol is not None:
#                kw = dict(ncol=legend_ncol)
#            else:
#                # TODO: get from legend label
#                legend_len = len(methods)
#                kw = dict(ncol=legend_len)
#            kw.update(dict(borderaxespad=0.0))
#            if title is not None:
#                kw.update(dict(title=title))
#            move_legend_facetgrid(g, **kw)
#    return g


# TODO: merge to set_move_legend()
def set_move_legend_ax(ax, legend):
    if legend is None:
        pass
    elif legend is False:
        remove_legend_sns(ax)
    else:
        move_legend_ax(ax)
    return ax


# TODO: somehow, legend of errorbar() do not have compatibility with sns.move_legend()
# create plt ver. ex. boosting_cpp_observe_wgt_snv_count::l.1409
def move_legend_ax(ax, **kwargs):
    # TODO: better way to add to dict
    if 'loc' not in kwargs:
        kwargs['loc'] = "upper left"
    if 'bbox_to_anchor' not in kwargs:
        kwargs['bbox_to_anchor'] = (1.05, 1)
    sns.move_legend(ax,
                    # should be in kwargs
                    # loc="upper left",
                    # bbox_to_anchor=(1.05, 1),
                    # title='Phenotype'
                    frameon=True,
                    fancybox=False,
                    edgecolor='black',
                    borderaxespad=0.0,
                    **kwargs,
                    )

    leg = get_legend_ax(ax)
    # change frame box line
    leg.get_frame().set_linewidth(1.0)
    return leg
