
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from logging import getLogger
from . import common as pcommon

# from ..validate.paper import config as pconfig
# from . import common as pcommon

logger = getLogger(__name__)

# text
color_text_gray_pptx = '#404040'


def plot_stat_heatmap(stats,
                      figtype=None,
                      # phe=None,
                      phed=None,
                      stat_name=None,
                      sex=None,
                      # palette=None,
                      statd=None,
                      ):

    xlabel = ylabel = 'Phenotype'

    stats=stats.copy()
    stats.index=stats.index.map(phed)
    stats.columns=stats.columns.map(phed)

    # figsize = (8, 8)
    _, figsize = pcommon.parse_figtype(figtype, True)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if stat_name == 'chisq':
        filt = (stats.values == 0.0)
        stats.values[filt] = np.nextafter(0, 1)
        logger.debug('stats: {}'.format(stats))

        annot_str = (np.asarray(['{:.2g}'.format(float(x)) if x >= 1e-300
                     else '<1e-300' for x in stats.to_numpy().flatten()])).reshape(stats.values.shape)

        # smaller, larger
        kwargs_heat = dict(
            cmap='Greens_r',
            #cmap='Blues_r',
                           vmin=0.0,
                           vmax=0.5,
                           norm=LogNorm(),
                           annot=annot_str,
                           fmt=''
                           )
    else:
        kwargs_heat = dict(
            cmap='Greens',
            #cmap='Blues',
                           annot=True)

    # only plot lower triangle
    mask = np.triu(np.ones_like(stats, dtype=bool))

    ax = sns.heatmap(data=stats,
                     # cmap='Blues',
                     # fmt='d',
                     cbar=False,
                     mask=mask,
                     ax=ax,
                     **kwargs_heat)

    ax.set(xlabel=xlabel, ylabel=ylabel)
    # ax.xaxis.set_label_position('top')

    ax.tick_params(left=False, bottom=False,
                   top=False,
                   #labelbottom=False,
                   #labeltop=True,
                   )


    ax.xaxis.set_tick_params(rotation=90)
    #ax.xaxis.set_tick_params(labelsize=15, rotation=90)

    title = statd[stat_name]
    if sex=='female':
        title+=' (female only)'
    ax.set_title(title)

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
