
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


def plot_scatter(stats,
                 figtype=None,
                 # phe=None,
                 phed=None,
                 stat_name=None,
                 ):

    _, figsize = pcommon.parse_figtype(figtype, True)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    ax = sns.scatterplot(stats, x='h2_liability', y='h2_liability_d',
                         )
    # ax = sns.pointplot(stats, x='h2_liability', y='h2_liability_d',
    #                   linestyles='none', errorbar=None)
    # linestyle='none', errorbar='none')

    # print(ax.get_xlim())
    # print(ax.get_ylim())
    lim_max = np.max((ax.get_xlim()[1], ax.get_ylim()[1]))
    lim_min = -0.005
    # lim_min = np.min((ax.get_xlim()[0], ax.get_ylim()[0]))
    # [print('lim', lim_max, lim_min)
    ax.set_xlim((lim_min, lim_max))
    ax.set_ylim((lim_min, lim_max))

    x = stats.loc[:, 'h2_liability']
    y = stats.loc[:, 'h2_liability_d']
    xerr = stats.loc[:, 'h2_liability_se'].to_numpy().reshape(1, -1)
    xerr = np.concatenate((xerr, xerr), axis=0)
    # [print('xerr')
    # [print(xerr.shape)
    yerr = stats.loc[:, 'h2_liability_se_d'].to_numpy().reshape(1, -1)
    yerr = np.concatenate((yerr, yerr), axis=0)
    # [print('yerr')
    # [print(yerr.shape)

    ax.errorbar(x, y, xerr=xerr, yerr=yerr,
                # 'none' does not plot (x,y) dot
                fmt='none',

                # vertical line width
                elinewidth=2,
                # ecolor=color
                )

    # annotation

    # y=x, y=x/10

    ax.set_xlabel('$h^2$')
    ax.set_ylabel('Dominance $h^2$')

    return fig
