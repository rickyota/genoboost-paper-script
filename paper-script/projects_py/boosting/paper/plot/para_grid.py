
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from logging import getLogger
from . import common as pcommon

# from ..validate.paper import config as pconfig
# from . import common as pcommon

logger = getLogger(__name__)

# text
color_text_gray_pptx = '#404040'


def plot_para_grid_heatmap(statd,
                           figtype=None,
                           method=None,
                           methodd=None,
                           phe=None,
                           phed=None,
                           stat_name=None,
                           # palette=None,
                           accd=None,
                           ):

    data = statd['grid']
    # row=#snv
    # col=learning rate

    # TODO: change col,index name of data
    # data.index = data.index.map(gmodeld)
    # data.columns = data.columns.map(gmodeld)

    xlabel = 'Learning rate'
    ylabel = '# SNVs'

    # figsize = (8, 8)
    _, figsize = pcommon.parse_figtype(figtype, True)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    ax = sns.heatmap(data=data,
                     cmap='Greens',
                     # cmap='Blues',
                     # fmt='d',
                     cbar=False,
                     annot=True,
                     ax=ax)

    # put frame on best para
    best_lr = statd['best_para']['lr']
    best_n = statd['best_para']['n']

    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.xaxis.set_label_position('top')

    ax.tick_params(left=False, bottom=False,
                   top=False,
                   labelbottom=False,
                   labeltop=True,
                   )
    for col_index, col_val in enumerate(data.columns.to_list()):
        for row_index, row_val in enumerate(data.index.to_list()):
            # float, int, np.float, np.int
            # print(type(col_val), type(row_val), type(best_lr), type(best_n))
            # print(col_val,row_val, best_lr,best_n)
            if col_val == best_lr and row_val == best_n:
                logger.debug('Put frame: {}, {}'.format(col_index, row_index))
                ax.add_patch(patches.Rectangle((col_index, row_index), 1, 1, fill=False, edgecolor='red', lw=3))
                # ax.add_patch(patches.Rectangle((col_index, row_index), 1, 1, fill=True, edgecolor='red', lw=10))

    title = accd[stat_name] + ': ' + methodd[method]
    # title = phed[phe] + ': ' + accd[stat_name]
    # TODO: use loc etc
    # or use fig.suptitle()
    ax.set_title(title)

    # fig.suptitle(title)
    # fig.suptitle(title,y=0.99)
    # pcommon.set_bottom_label(fig, title, 0.02)

    # change axis label color to gray
    # ax.spines[['left', 'bottom']].set_color(color_text_gray_pptx)
    # ax.xaxis.label.set_color(color_text_gray_pptx)
    # ax.yaxis.label.set_color(color_text_gray_pptx)
    # ax.tick_params(axis='x', colors=color_text_gray_pptx)
    # ax.tick_params(axis='y', colors=color_text_gray_pptx)

    # change tick label color
    # palettep = [palette[k] for k in gmodels]
    # [t.set_color(c) for t, c in zip(ax.xaxis.get_ticklabels(), palettep)]
    # [t.set_color(c) for t, c in zip(ax.yaxis.get_ticklabels(), palettep)]

    fig.tight_layout()

    return fig


def plot_para_grid_lineplot(statd,
                            figtype=None,
                            method=None,
                            methodd=None,
                            phe=None,
                            phed=None,
                            stat_name=None,
                            # palette=None,
                            accd=None,
                            ):

    data = statd['allstat']
    # data = statd['grid']
    # row=#snv
    # col=learning rate

    # TODO: change col,index name of data
    # data.index = data.index.map(gmodeld)
    # data.columns = data.columns.map(gmodeld)

    # xlabel = 'Learning rate'
    xlabel = '# SNVs'
    ylabel = accd[stat_name]

    # figsize = (8, 8)
    _, figsize = pcommon.parse_figtype(figtype, True)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    logger.debug('data: {}'.format(data))

    ax = sns.lineplot(data=data,
                      x='nsnv_use',
                      y=stat_name,
                      hue='lr',
                      style='lr',
                      palette='coolwarm',
                      # palette='cool',
                      ax=ax,
                      markers='o',
                      markersize=8,
                      #markersize=10,
                      )

    # put marker on best para
    best_lr = statd['best_para']['lr']
    best_n = statd['best_para']['n']
    best_stat = data.loc[(data['lr'] == best_lr) & (data['n'] == best_n), stat_name]
    assert len(best_stat) == 1

    ax.plot(best_n, best_stat, 'X', color='k', markersize=20, markeredgewidth=0.0)
    # linewidth no effect
    # ax.plot(best_n, best_stat, 'X', color='k', markersize=20, linewidth=1.5)
    # ax.plot(best_n, best_stat, 'ko', markersize=20, fillstyle='none', markeredgewidth=1.5)

    # ax = sns.heatmap(data=data,
    #                 cmap='Greens',
    #                 # cmap='Blues',
    #                 # fmt='d',
    #                 cbar=False,
    #                 annot=True,
    #                 ax=ax)

    ax.set(xlabel=xlabel, ylabel=ylabel)
    # ax.xaxis.set_label_position('top')

    # ax.tick_params(left=False, bottom=False,
    #               top=False,
    #               labelbottom=False,
    #               labeltop=True,
    #               )

    title = methodd[method]
    # title = accd[stat_name] + ': ' + methodd[method]
    # title = phed[phe] + ': ' + accd[stat_name]
    # TODO: use loc etc
    # or use fig.suptitle()
    ax.set_title(title)

    # ax.legend(loc='lower left', bbox_to_anchor=(1., 0))
    pcommon.move_legend_ax(
        ax, loc='lower left', bbox_to_anchor=(1.05, 0),
        title='learning rate')

    # fig.suptitle(title)
    # fig.suptitle(title,y=0.99)
    # pcommon.set_bottom_label(fig, title, 0.02)

    # change axis label color to gray
    # ax.spines[['left', 'bottom']].set_color(color_text_gray_pptx)
    # ax.xaxis.label.set_color(color_text_gray_pptx)
    # ax.yaxis.label.set_color(color_text_gray_pptx)
    # ax.tick_params(axis='x', colors=color_text_gray_pptx)
    # ax.tick_params(axis='y', colors=color_text_gray_pptx)

    # change tick label color
    # palettep = [palette[k] for k in gmodels]
    # [t.set_color(c) for t, c in zip(ax.xaxis.get_ticklabels(), palettep)]
    # [t.set_color(c) for t, c in zip(ax.yaxis.get_ticklabels(), palettep)]

    fig.tight_layout()

    return fig
