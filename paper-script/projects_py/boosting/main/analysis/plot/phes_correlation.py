
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from logging import getLogger


logger = getLogger(__name__)


def plot_phe_correlation(droot, stats, stat_name):
    figsize = (12, 12)
    # figsize = (8, 8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if stat_name == 'chisq':
        filt = (stats.values == 0.0)
        stats.values[filt] = np.nextafter(0, 1)
        print('stats', stats)
        # smaller, larger

        # annot_str = (np.asarray(['{:.2g}'.format(float(x)) if x != np.nextafter(0, 1)
        annot_str = (np.asarray(['{:.2g}'.format(float(x)) if x >= 1e-300
                     else '<1e-300' for x in stats.to_numpy().flatten()])).reshape(stats.values.shape)
        # print(annot_str)

        kwargs_heat = dict(cmap='Blues_r',
                           vmin=0.0,
                           vmax=0.5,
                           norm=LogNorm(),
                           annot=annot_str,
                           fmt=''
                           )
    else:
        kwargs_heat = dict(cmap='Blues',
                           annot=True)

    mask = np.triu(np.ones_like(stats, dtype=bool))
    # mask=~np.tri(len(stats),k=-1,dtype=bool)

    ax = sns.heatmap(data=stats,
                     # cmap=cmap,
                     # cmap='Blues',
                     # fmt='d',
                     cbar=False,
                     # annot=True,
                     mask=mask,
                     ax=ax,
                     **kwargs_heat)

    # if stat_name == 'chisq':
    #    # adjust 0.0
    #    pass
    #    for i in range(len(stats)):
    #        for j in range(len(stats)):
    #            stat = stats.iloc[i, j]
    #            if stat == 0.0:
    #                #stat = np.nextafter(0, 1)
    #                stat_str='<1e-300'
    #            else:
    #                stat_str='.2g'.format(stat)
    #            _ = ax.text(j, i, stat_str)

    ax.set_title(stat_name)

    fig.tight_layout()
    fname = droot / (stat_name + '.heatmap.figtype' + '0')
    logger.debug('fname {}'.format(fname))
    fig.savefig(fname.parent / (fname.name + '.png'))
    plt.close()
