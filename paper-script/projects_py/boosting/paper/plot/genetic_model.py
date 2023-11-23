
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from logging import getLogger

#from ..validate.paper import config as pconfig
from .. import config as pconfig
from . import common as pcommon

logger = getLogger(__name__)

# text
color_text_gray_pptx = '#404040'


def plot_genetic_model_gwas_heatmap_gwasmethods(ss_wgt, gwas_model_est='gwas'):
    # row=GenoBoost est model
    # col=GWAS est model
    # y=GB effect of the model
    # x=GWAS effect of model
    # ** only plot the estimated plot in the row **

    # cannot use relplot since x and y are different

    if gwas_model_est == 'gwas':
        gwas_model_est = 'model_est_gwas'
        gwas_model_est_plot = gwas_model_est.capitalize()
        # gwas_model_est_plot = 'Model_est_gwas'
        xlabel = 'GWAS'
    elif gwas_model_est == 'aic':
        gwas_model_est = 'model_est_aic'
        gwas_model_est_plot = gwas_model_est.capitalize()
        # gwas_model_est_plot = 'Model_est_aic'
        xlabel = 'AIC'
    else:
        raise RuntimeError

    ss_wgt = ss_wgt.copy()
    ss_wgt.loc[:, 'Boost_model_est'] = ss_wgt['boost_model_est'].str.title()
    ss_wgt.loc[:, gwas_model_est_plot] = ss_wgt[gwas_model_est].str.title()
    # ss_wgt.loc[:,'model_est_gwas']= ss_wgt['model_est_gwas'].str.title()

    gmodels = ['Overdom', 'Dom', 'Add', 'Rec', 'Overrec']
    # gmodels = ['overdom', 'dom', 'add', 'rec', 'overrec']

    figsize = (8, 8)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    data = ss_wgt[['Boost_model_est', gwas_model_est_plot]].copy()
    # data = ss_wgt[['boost_model_est', 'model_est_gwas']].copy()
    data.loc[:, 'dummy'] = 1
    print('data', data.head())

    data = pd.pivot_table(data,
                          values='dummy',
                          index='Boost_model_est', columns=gwas_model_est_plot,
                          # index='boost_model_est', columns='model_est_gwas',
                          aggfunc='count')
    print('data pivot', data)

    def fill_col_index(data, gmodels):
        for gmodel in gmodels:
            # fill index
            if gmodel not in data.index:
                data.loc[gmodel, :] = np.nan
            # fill col
            if gmodel not in data:
                data.loc[:, gmodel] = np.nan
        return data
    data = fill_col_index(data, gmodels)

    data = data.loc[gmodels, gmodels]

    # FIXME: use categorical type for boost_model_est
    # then, these are not necessary
    data = data.fillna(0.0)
    data = data.astype(int)
    print('data pivot', data)

    ax = sns.heatmap(data=data,
                     cmap='Blues',
                     fmt='d',
                     cbar=False,
                     annot=True,
                     ax=ax)

    # ax.set(xlabel='GWAS',
    ax.set(xlabel=xlabel, ylabel='GenoBoost')

    ax.tick_params(left=False, bottom=False)

    ax.spines[['left', 'bottom']].set_color(color_text_gray_pptx)
    ax.xaxis.label.set_color(color_text_gray_pptx)
    ax.yaxis.label.set_color(color_text_gray_pptx)
    ax.tick_params(axis='x', colors=color_text_gray_pptx)
    ax.tick_params(axis='y', colors=color_text_gray_pptx)

    cb = sns.color_palette("colorblind", 10)
    # ['overdom', 'dom','add','rec', 'overrec']
    # should be different colors as odds ratio
    mypalette4 = [
        cb[6], cb[4], cb[2], cb[0], cb[9]
    ]

    [t.set_color(c) for t, c in zip(ax.xaxis.get_ticklabels(), mypalette4)]
    [t.set_color(c) for t, c in zip(ax.yaxis.get_ticklabels(), mypalette4)]

    fig.tight_layout()

    return fig


def plot_genetic_model_gwas_heatmap(data, gwas_model_est_criteria='gwas',
                                    # gmodels=None,
                                    gmodeld=None,
                                    palette=None,
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

    if gwas_model_est_criteria == 'gwas':
        # gwas_model_est = 'model_est_gwas'
        # gwas_model_est_plot = 'Model_est_gwas'
        xlabel = 'GWAS'
    elif gwas_model_est_criteria == 'aic':
        # gwas_model_est = 'model_est_aic'
        # gwas_model_est_plot = 'Model_est_aic'
        xlabel = 'AIC'

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

    ax.set(xlabel=xlabel, ylabel='GenoBoost')

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


def genetic_model_prop(
        stats,
        figtype=None,
        xlabel_fontsize=None,
        ylabel_fontsize=None,
        xtick_fontsize=None,
        ytick_fontsize=None,
):
    _, fig_aspect = pcommon.parse_figtype(figtype)
    figsize = (fig_aspect['height'], fig_aspect['height'] * fig_aspect['aspect'])

    """
        # hatch on barplot
        https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
        https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
        """

    genetic_models = stats

    # ['overdom', 'dom','add','rec', 'overrec']
    # should be different colors as odds ratio
    cb = sns.color_palette("colorblind", 10)
    mypalette4 = [
        cb[6], cb[4], cb[2], cb[0], cb[9]
    ]
    # for stacked bar
    mypalette4r = list(reversed(mypalette4))

    # first, bottom
    # change legend order later
    columns_order = ['overrec', 'rec', 'add', 'dom', 'overdom']
    # columns_order = ['overdom', 'dom', 'add', 'rec', 'overrec']
    # columns_order = ['add', 'overdom', 'dom', 'overrec', 'rec']
    genetic_models = genetic_models.loc[:, columns_order]
    genetic_models.index = genetic_models.index.map(pconfig._phed_capital3())
    # genetic_models.index = genetic_models.index.map(get_phed())
    genetic_models.columns = genetic_models.columns.map(pconfig.genetic_modeld())
    hatchs = ['//', '\\\\', '', '||', '--', ]
    # hatchs = ['//', '\\\\', '..', '||', '--', ]
    # hatchs = ['--', '++', '.', '//', '\\']
    # order is inversed
    mypalette3_alpha = [v + (1.0,) for v in mypalette4r]
    color_hatch = dict(zip(mypalette3_alpha, hatchs))
    with sns.plotting_context('poster'):
        with sns.color_palette(mypalette4r, 5):
            print('gm', genetic_models)
            genetic_models_percent = genetic_models.copy()
            for index in genetic_models_percent.index:
                print('index', index)
                print('row', genetic_models_percent.loc[index, :])
                genetic_models_percent.loc[index, :] = genetic_models_percent.loc[index, :] * 100 / genetic_models_percent.loc[index, :].sum()
            print('genetic_models_percent', genetic_models_percent)

            fig, ax = plt.subplots(1, 1, figsize=figsize)
            # fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            # edgecolor is for hatch color but also change edge
            ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False, edgecolor='w')
            # ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False)
            for mi, bar in enumerate(ax.patches):
                hatch = color_hatch[bar.get_facecolor()]
                bar.set_hatch(hatch)
            ax.set_xlabel('Phenotype', fontsize=xlabel_fontsize)
            ax.set_ylabel('Proportion of SNVs', fontsize=ylabel_fontsize)
            # ax.set(xlabel='Phenotype',
            # ylabel='Proportion of SNVs')
            # ax.set_xlabel('Phenotype', fontsize=28)
            # ax.set_ylabel('Proportion of SNVs', fontsize=28)

            ax.set_ylim((0, 100))
            ax.set_yticks((0, 25, 50, 75, 100))
            ax.set_yticklabels(('0%', '25%', '50%', '75%', '100%'))
            # https://analytics-note.xyz/programming/matplotlib-tick-rotation/
            if xtick_fontsize is None:
                xtick_fontsize = 20
            ax.xaxis.set_tick_params(labelsize=xtick_fontsize, rotation=90)
            # set_tick_params() does not have linespacing
            plt.xticks(linespacing=1.0)


            if ytick_fontsize is None:
                ytick_fontsize = 15
            ax.yaxis.set_tick_params(labelsize=ytick_fontsize)

            # fig.tight_layout()
            # fig.savefig(fout + '.percent.bar.nolegend.xrot.png')
            # fig.savefig(fout + '.percent.bar.nolegend.xrot.pdf')

            return fig
