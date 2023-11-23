
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from logging import getLogger

# from ..validate.paper import config as pconfig
# from . import common as pcommon

logger = getLogger(__name__)

# text
color_text_gray_pptx = '#404040'


def plot_genetic_model_gwas_heatmap(data,
                                    xy,
                                    # gmodels=None,
                                    gmodeld=None,
                                    palette=None,
                                    title=None,
                                    ):
    # row=GenoBoost est model
    # col=GWAS est model
    # y=GB effect of the model
    # x=GWAS effect of model
    # ** only plot the estimated plot in the row **

    # cannot use relplot since x and y are different

    # if gmodels is None:
    # gmodels = ['overdom', 'dom', 'add', 'rec', 'overrec']

    if data.index.tolist() != data.columns.tolist():
        logger.error('index vs cols: {} vs {}'.format(data.index, data.columns))
        raise RuntimeError

    gmodels = list(data.index)

    if xy == ('gwas', 'boost'):
        xlabel, ylabel = 'GWAS', 'GenoBoost'
    elif xy == ('true', 'boost'):
        xlabel, ylabel = 'True model', 'GenoBoost'
    elif xy == ('true', 'gwas'):
        xlabel, ylabel = 'True model', 'GWAS'
    else:
        raise RuntimeError

    # TODO: change col,index name of data
    data.index = data.index.map(gmodeld)
    data.columns = data.columns.map(gmodeld)

    figsize = (8, 8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    ax = sns.heatmap(data=data,
                     cmap='Blues',
                     fmt='d',
                     cbar=False,
                     annot=True,
                     ax=ax)

    ax.set_title(title)
    ax.set(xlabel=xlabel, ylabel=ylabel)

    ax.tick_params(left=False, bottom=False)

    # change axis label color to gray
    ax.spines[['left', 'bottom']].set_color(color_text_gray_pptx)
    ax.xaxis.label.set_color(color_text_gray_pptx)
    ax.yaxis.label.set_color(color_text_gray_pptx)
    ax.tick_params(axis='x', colors=color_text_gray_pptx)
    ax.tick_params(axis='y', colors=color_text_gray_pptx)

    # change tick label color
    palettep = [palette[k] for k in gmodels]
    [t.set_color(c) for t, c in zip(ax.xaxis.get_ticklabels(), palettep)]
    [t.set_color(c) for t, c in zip(ax.yaxis.get_ticklabels(), palettep)]

    fig.tight_layout()

    return fig


def plot_genetic_model(
    genetic_models,
    gmodeld,
    title,
):
    """
    # hatch on barplot
    https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
    https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
    """

    phed = {
        'true_model_est': 'True model',
        'model_est_gwas': 'GWAS',
        'boost_model_est': 'GenoBoost'
    }

    cb = sns.color_palette("colorblind", 10)
    #  ['overdom', 'dom','add','rec', 'overrec']
    #  should be different colors from odds ratio
    mypalette4 = [
        cb[6], cb[4], cb[2], cb[0], cb[9]
    ]
    #  for stacked bar
    mypalette4r = list(reversed(mypalette4))

    # first, bottom
    # change legend order later
    columns_order = ['overrec', 'rec', 'add', 'dom', 'overdom']
    # columns_order = ['overdom', 'dom', 'add', 'rec', 'overrec']
    # columns_order = ['add', 'overdom', 'dom', 'overrec', 'rec']
    genetic_models = genetic_models.loc[:, columns_order]
    genetic_models.index = genetic_models.index.map(phed)
    # genetic_models.index = genetic_models.index.map(get_phed())
    genetic_models.columns = genetic_models.columns.map(gmodeld)

    hatchs = ['//', '\\\\', '', '||', '--', ]
    # hatchs = ['//', '\\\\', '..', '||', '--', ]
    # hatchs = ['--', '++', '.', '//', '\\']
    # order is inversed
    mypalette3_alpha = [v + (1.0,) for v in mypalette4r]
    color_hatch = dict(zip(mypalette3_alpha, hatchs))

    # with sns.plotting_context('poster'):
    with sns.color_palette(mypalette4r, 5):
        logger.debug('gm {}'.format(genetic_models))
        genetic_models_percent = genetic_models.copy()
        for index in genetic_models_percent.index:
            # ng: when success, this is series; distinguish by series or dataframe
            # if len(genetic_models_percent.loc[index, :].index) > 1:
            #    logger.error('index: {}'.format(index))
            #    logger.error('stat {}'.format(genetic_models_percent.loc[index, :].index))
            #    logger.error('len {}'.format(len(genetic_models_percent.loc[index, :].index)))
            #    # index is nan
            #    raise RuntimeError
            # could be string
            # if np.isnan(index):
            #    logger.error('index: {}'.format(index))
            #    raise RuntimeError('index is NaN')
            logger.debug('index: {}'.format(index))
            logger.debug('row: {}'.format(genetic_models_percent.loc[index, :]))
            # print('index', index)
            # print('row', genetic_models_percent.loc[index, :])
            genetic_models_percent.loc[index, :] = 100 * genetic_models_percent.loc[index, :] / genetic_models_percent.loc[index, :].sum()
        logger.debug('genetic_models_percent: {}'.format(genetic_models_percent))

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        # edgecolor is for hatch color but also change edge
        ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False, edgecolor='w')
        # ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False)
        for mi, bar in enumerate(ax.patches):
            hatch = color_hatch[bar.get_facecolor()]
            bar.set_hatch(hatch)

        ax.set_title(title, fontsize=15)
        ax.set(xlabel='',
               ylabel='Proportion of SNVs')

        ax.set_ylim((0, 100))
        ax.set_yticks((0, 25, 50, 75, 100))
        ax.set_yticklabels(('0%', '25%', '50%', '75%', '100%'))
        # https://analytics-note.xyz/programming/matplotlib-tick-rotation/
        ax.yaxis.set_tick_params(labelsize=20)
        ax.xaxis.set_tick_params(labelsize=15)
        # ax.xaxis.set_tick_params(labelsize=15, rotation=90)
        fig.tight_layout()
        # fig.savefig(fout + '.percent.bar.nolegend.xrot.png')
        # fig.savefig(fout + '.percent.bar.nolegend.xrot.pdf')
        return fig


# def validate_genetic_model_old(droot, pformats, stats, phes):
#    # TMP cvi=0 only
#    gmodels = ['overrec', 'rec', 'add', 'dom', 'overdom']

#    def stats_gmodels(stats, cvi):
#        logger.debug('stats {}'.format(stats))
#        logger.debug('stats {}'.format(stats.info()))
#        stats_cvi = stats.loc[stats['cvi'] == cvi, ['phe'] + gmodels]
#        # check if only one stat for each phe
#        assert (~stats_cvi['phe'].duplicated()).all()
#        gmodel = stats_cvi.set_index('phe')
#        return gmodel

#    cvi = 0
#    genetic_models = stats_gmodels(stats, cvi)
#    # genetic_models = stats

#    phes_ext = pconfig.phes_12()
#    genetic_models = genetic_models.loc[phes_ext, :]

#    # all
#    if 'paper' in pformats:
#        """
#        # hatch on barplot
#        https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
#        https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
#        """

#        cb = sns.color_palette("colorblind", 10)
#        # ['overdom', 'dom','add','rec', 'overrec']
#        # should be different colors as odds ratio
#        mypalette4 = [
#            cb[6], cb[4], cb[2], cb[0], cb[9]
#        ]
#        # for stacked bar
#        mypalette4r = list(reversed(mypalette4))

#        # first, bottom
#        # change legend order later
#        columns_order = ['overrec', 'rec', 'add', 'dom', 'overdom']
#        # columns_order = ['overdom', 'dom', 'add', 'rec', 'overrec']
#        # columns_order = ['add', 'overdom', 'dom', 'overrec', 'rec']
#        genetic_models = genetic_models.loc[:, columns_order]
#        genetic_models.index = genetic_models.index.map(pconfig.phed_capital3())
#        # genetic_models.index = genetic_models.index.map(get_phed())
#        genetic_models.columns = genetic_models.columns.map(pconfig.genetic_modeld())
#        hatchs = ['//', '\\\\', '', '||', '--', ]
#        # hatchs = ['//', '\\\\', '..', '||', '--', ]
#        # hatchs = ['--', '++', '.', '//', '\\']
#        # order is inversed
#        mypalette3_alpha = [v + (1.0,) for v in mypalette4r]
#        color_hatch = dict(zip(mypalette3_alpha, hatchs))
#        with sns.plotting_context('poster'):
#            with sns.color_palette(mypalette4r, 5):
#                print('gm', genetic_models)
#                genetic_models_percent = genetic_models.copy()
#                for index in genetic_models_percent.index:
#                    print('index', index)
#                    print('row', genetic_models_percent.loc[index, :])
#                    genetic_models_percent.loc[index, :] = 100 * genetic_models_percent.loc[index, :] / genetic_models_percent.loc[index, :].sum()
#                print('genetic_models_percent', genetic_models_percent)

#                fig, ax = plt.subplots(1, 1, figsize=(8, 8))
#                # edgecolor is for hatch color but also change edge
#                ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False, edgecolor='w')
#                # ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False)
#                for mi, bar in enumerate(ax.patches):
#                    hatch = color_hatch[bar.get_facecolor()]
#                    bar.set_hatch(hatch)
#                ax.set(xlabel='Phenotype',
#                       ylabel='Proportion of SNVs')

#                ax.set_ylim((0, 100))
#                ax.set_yticks((0, 25, 50, 75, 100))
#                ax.set_yticklabels(('0%', '25%', '50%', '75%', '100%'))
#                # https://analytics-note.xyz/programming/matplotlib-tick-rotation/
#                ax.yaxis.set_tick_params(labelsize=20)
#                ax.xaxis.set_tick_params(labelsize=15, rotation=90)
#                fig.tight_layout()
#                # fig.savefig(fout + '.percent.bar.nolegend.xrot.png')
#                # fig.savefig(fout + '.percent.bar.nolegend.xrot.pdf')

#                figtype = '1-8-1'

#                fname = droot / ('genetic-model.12phe.percent_bar.nolegend.paper.gwasmethods.figtype' + figtype)
#                pio.savefig_png_svg(fig, fname)
#                plt.close()
