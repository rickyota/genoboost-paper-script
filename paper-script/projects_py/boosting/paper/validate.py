# [seaborn objects](https://seaborn.pydata.org/tutorial/objects_interface.html)
# Should I split validate/paper/auc.py and plot/auc.py

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
# from collections import defaultdict
from sklearn.linear_model import LinearRegression
import math
from logging import getLogger

from .plot import acc as plot_acc
from .plot import odds_ratio as plot_odds_ratio
from .plot import nsnv as plot_nsnv
from .plot import genetic_model as plot_genetic_model
from .plot import para_grid as plot_para_grid
from .plot import phes_correlation as plot_phes_correlation
from .plot import manhattan as plot_manhattan
from .plot import sample_score as plot_sample_score
from .plot import h2 as plot_h2
from .plot import io as pio
from .plot import config as pltconfig
# from .plot import common as pcommon
from . import config as pconfig
from . import public_stat


logger = getLogger(__name__)

# or local?
# but even for local, once called, the latter codes reflect it
pltconfig.plt_font_paper()


def extract_boosing_methods(methods, ext):
    # stay the order in methods same
    methods_ext = []
    for method in methods:
        if method.startswith('boosting_') and method not in ext:
            continue
        else:
            methods_ext.append(method)
    return methods_ext


# TODO: stat_name => plot_on
def validate_para_grid(droot, pformats, stat, stat_name, methodm, phe):

    # stat_heatmap = create_heatmap_table(ss_wgt, gwas_model_est_criteria, gmodels)

    if 'paper' in pformats:
        with sns.plotting_context('poster',  # font_scale=1.0,
                                  rc={'font.size': 10,
                                      'axes.titlesize': 22, 'axes.labelsize': 18,
                                      'ytick.labelsize': 16, 'xtick.labelsize': 16}
                                  ):
            # figtype = '0-12-0.75'
            figtype = '0-12-0.6'
            fig = plot_para_grid.plot_para_grid_heatmap(stat,
                                                        figtype=figtype,
                                                        method=methodm,
                                                        methodd=pconfig.methodd([methodm], 'capital'),
                                                        phe=phe,
                                                        phed=pconfig.phed([phe], 'capital'),
                                                        # palette=pconfig.palette_hls_2_genetic_model(),
                                                        stat_name=stat_name,
                                                        accd=pconfig.accd(),
                                                        )
            fig.tight_layout()
            # fig.tight_layout(h_pad=-0.8, w_pad=0)
            fname = droot / ('heatmap' + '.paper.figtype' + figtype)
            logger.debug('fname {}'.format(fname))
            pio.savefig_png_svg(fig, fname, dpi=1200)
            plt.close()

    if 'paper' in pformats:
        with sns.plotting_context('poster',  # font_scale=1.0,
                                  rc={'font.size': 10,
                                      'axes.titlesize': 22, 'axes.labelsize': 18,
                                      'ytick.labelsize': 16, 'xtick.labelsize': 16}
                                  ):
            # figtype = '0-12-0.75'
            # figtype = '0-12-0.6'
            figtype = '0-6-2'
            fig = plot_para_grid.plot_para_grid_lineplot(stat,
                                                         figtype=figtype,
                                                         method=methodm,
                                                         methodd=pconfig.methodd([methodm], 'capital'),
                                                         phe=phe,
                                                         phed=pconfig.phed([phe], 'capital'),
                                                         # palette=pconfig.palette_hls_2_genetic_model(),
                                                         stat_name=stat_name,
                                                         accd=pconfig.accd(),
                                                         )
            fig.tight_layout()
            # fig.tight_layout(h_pad=-0.8, w_pad=0)
            fname = droot / ('lineplot' + '.paper.figtype' + figtype)
            logger.debug('fname {}'.format(fname))
            pio.savefig_png_svg(fig, fname, dpi=1200)
            plt.close()


# file name must include plot_on
def validate_acc(droot, pformats, stats, plot_on, methods, methodd_type, phes):

    # all
    if 'paper' in pformats:
        # width=180mm=7inch
        # 5inch(default) * 6 =30inch
        # scale=7/30=1/4.3
        #
        # minimum font ytick.labelsize=18
        # 18pt*1/4.3=4.2pt
        # shold make this size >5pt as instructed by Nature
        # -> but looks fine in pptx
        #

        # if len(methods) >= 10:
        #    figtypes = ['3-5-0.8', '4-5-1.5']
        # else:
        #    figtypes = ['3-5-0.8']
        # figtypes = ['3-5-0.8', '4-5-1.5']
        figtypes = ['3-5-0.8', '4-5-1.5', '5-5-0.9', '6-5-1']

        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                    'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                      methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
                                                      methodd=pconfig.methodd(methods, 'capital2_roman'),
                                                      phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                      col_wrap=min(6, len(phes)),
                                                      palette=pconfig.palette_hls_2_methods(),
                                                      sharey=None,
                                                      ylabel=pconfig.accd()[plot_on],
                                                      xtick_label_format='norot',
                                                      # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                      xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
                                                      ytick_ntick=4,
                                                      legend='horizontal')
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / (plot_on + '.all.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

                # vlines to split boosting and others
                fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                      methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
                                                      methodd=pconfig.methodd(methods, 'capital2_roman'),
                                                      phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                      col_wrap=min(6, len(phes)),
                                                      palette=pconfig.palette_hls_2_methods(),
                                                      sharey=None,
                                                      ylabel=pconfig.accd()[plot_on],
                                                      xtick_label_format='norot',
                                                      # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                      xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
                                                      ytick_ntick=4,
                                                      legend='horizontal',
                                                      vlines=[2.5],
                                                      )
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / (plot_on + '.all.vline.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

        # ex-chr6, noneur
        # - align legend of snpnet and snpnet_noneur
        # - add middle vline
        # if methodd_type exists
        if len(methodd_type.values()) > 0:
            # col_wrap
            figtypes = ['4-5-1.5', '5-5-1.7', '6-5-1.8']
            for figtype in figtypes:
                with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                        'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                    fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                          methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
                                                          methodd=pconfig.methodd(methods, 'capital2_roman'),
                                                          phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                          col_wrap=min(4, len(phes)),
                                                          palette=pconfig.palette_hls_2_methods(),
                                                          sharey=None,
                                                          ylabel=pconfig.accd()[plot_on],
                                                          xtick_label_format='norot',
                                                          # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                          xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
                                                          ytick_ntick=4,
                                                          legend='horizontal',
                                                          legend_ncol=math.ceil(len(methods) / 2),
                                                          vlines=[len(methods) / 2 - 0.5],
                                                          )
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.all.col_wrap-4.paper.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

                    if len(phes) == 5:
                        # col_wrap
                        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                            fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                                  methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
                                                                  methodd=pconfig.methodd(methods, 'capital2_roman'),
                                                                  phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                                  col_wrap=3,
                                                                  palette=pconfig.palette_hls_2_methods(),
                                                                  sharey=None,
                                                                  ylabel=pconfig.accd()[plot_on],
                                                                  xtick_label_format='norot',
                                                                  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                  xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
                                                                  ytick_ntick=4,
                                                                  legend='horizontal',
                                                                  legend_ncol=math.ceil(len(methods) / 2),
                                                                  vlines=[len(methods) / 2 - 0.5],
                                                                  )
                            # 0.8 time fontheight
                            fig.tight_layout(h_pad=0.8, w_pad=0)

                            fname = droot / (plot_on + '.all.col_wrap-3.paper.figtype' + figtype)
                            logger.debug('fplot: {}'.format(fname))
                            pio.savefig_png_svg(fig, fname)
                            plt.close()

    # 12phe should be in script
    # integ, add, nonadd + previous methods, 12phe
    # if 'paper' in pformats:
    #    figtype = '3-5-0.8'
    #    with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
    #                                            'xtick.major.width': 1.5, 'xtick.major.size': 8}):
    #        # if {'boosting_integrate', 'boosting_add', 'boosting_nonadd'} in set(methods):
    #        methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_nonadd', 'boosting_add'])
    #        fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
    #                                              methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
    #                                              phes=pconfig.phes_12(), phed=pconfig.phed(pconfig.phes_12(), 'capital2'),
    #                                              col_wrap=6,
    #                                              palette=pconfig.palette_hls_2_methods(),
    #                                              sharey=None,
    #                                              ylabel=pconfig.accd()[plot_on],
    #                                              xtick='norot',
    #                                              xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
    #                                              ytick_ntick=4,
    #                                              legend='horizontal',
    #                                              vlines=[2.5])
    #        # 0.8 time fontheight
    #        fig.tight_layout(h_pad=0.8, w_pad=0)

    #        fname = droot / (plot_on + '.integ-add-nonadd-prev_12phe.paper.figtype' + figtype)
    #        logger.debug('fplot: {}'.format(fname))
    #        pio.savefig_png_svg(fig, fname)
    #        plt.close()

    # integ + prev
    if 'paper' in pformats:
        # figtype = '8-5-0.7'
        figtype = '9-5-0.8'
        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            # if {'boosting_integrate', 'boosting_add', 'boosting_nonadd'} <= set(methods):
            if {'boosting_integrate'} <= set(methods):
                methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
                fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                      methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                      phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                      col_wrap=6,
                                                      palette=pconfig.palette_hls_2_methods(),
                                                      sharey=None,
                                                      ylabel=pconfig.accd()[plot_on],
                                                      xtick_label_format='norot',
                                                      xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                      ytick_ntick=4,
                                                      legend='horizontal')
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / (plot_on + '.integ-prev.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

    if 'poster' in pformats:
        # large 4 phes and small 5 phes (atm, ad, bc, cad t2d)
        # large 4 phes
        # figtypes = ['1-5-0.8', '2-4-0.8']
        figtypes = ['2-4-0.8', '3-5-0.6']
        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                    'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                # if {'boosting_integrate', 'boosting_add', 'boosting_nonadd'} <= set(methods):
                if {'boosting_integrate'} <= set(methods):
                    methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
                    phes_ext = ['ra', 'psr', 'gout', 'ibd']
                    fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                          methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital'),
                                                          phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital3'),
                                                          col_wrap=len(phes_ext),
                                                          palette=pconfig.palette_hls_2_methods(),
                                                          sharey=None,
                                                          ylabel=pconfig.accd()[plot_on],
                                                          xtick_format='none',
                                                          xtick_label_format='none',
                                                          ytick_ntick=3,
                                                          legend='horizontal')
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-prev_phes4.poster.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

        # large 5 phes
        # figtypes = ['1-3.5-1']
        figtypes = ['2-4-0.8']
        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                    'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                # if {'boosting_integrate', 'boosting_add', 'boosting_nonadd'} <= set(methods):
                if {'boosting_integrate'} <= set(methods):
                    methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
                    phes_ext = ['atm', 'ad', 'bc', 'cad', 't2d']
                    fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                          methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital'),
                                                          phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital3'),
                                                          col_wrap=len(phes_ext),
                                                          palette=pconfig.palette_hls_2_methods(),
                                                          sharey=None,
                                                          ylabel=pconfig.accd()[plot_on],
                                                          xtick_format='none',
                                                          xtick_label_format='none',
                                                          ytick_ntick=3,
                                                          legend='horizontal')
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-prev_phes5.poster.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

    # integ + prev; modify xtick
    if 'paper' in pformats:
        # figtype = '8-5-0.7'
        figtype = '9-5-0.8'
        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            # if {'boosting_integrate', 'boosting_add', 'boosting_nonadd'} <= set(methods):
            if {'boosting_integrate'} <= set(methods):
                methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
                xtick_label = [x if x != 'viii' else '  viii' for x in pconfig.create_xtick_labels(methods_ext, 'roman_number')]
                fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                      methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                      phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                      col_wrap=6,
                                                      palette=pconfig.palette_hls_2_methods(),
                                                      sharey=None,
                                                      ylabel=pconfig.accd()[plot_on],
                                                      xtick_label_format='norot',
                                                      xtick_label=xtick_label,
                                                      ytick_ntick=4,
                                                      legend='horizontal')
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / (plot_on + '.integ-prev.xtick.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

    # x TODO: create .sh for phes2 -> but this requires a lot of sh files
    # extracting methods can be here
    # integ, add, nonadd only, phes2=['ra', 'atm', 'psr']
    if 'paper' in pformats:
        figtype = '5-5-1.0'
        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            methods_ext = ['boosting_integrate', 'boosting_nonadd', 'boosting_add']
            if set(methods_ext) <= set(methods):
                # phes_ext = ['ra', 'atm', 'psr']
                phes_ext = ['ra', 'psr', 'atm']
                fig = plot_acc.plot_acc_box_strip_cat(
                    stats,
                    plot_on,
                    figtype=figtype,
                    methods=methods_ext,
                    methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                    phes=phes_ext,
                    phed=pconfig.phed(phes_ext, 'capital2'),
                    col_wrap=1,
                    palette=pconfig.palette_hls_2_methods(),
                    sharey=None,
                    ylabel=pconfig.accd()[plot_on],
                    xtick_label_format='norot',
                    xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                    ytick_ntick=4,
                    legend='horizontal',
                    legend_ncol=1)
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / (plot_on + '.integ-add-nonadd_phes2.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

    if 'paper' in pformats:
        figtype = '5-5-1.0'
        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            methods_ext = ['boosting_integrate', 'boosting_nonadd', 'boosting_add']
            if set(methods_ext) <= set(methods):
                # phes_ext = ['ra', 'atm', 'psr']
                phes_ext = ['ra', 'psr', 'gout']
                fig = plot_acc.plot_acc_box_strip_cat(
                    stats,
                    plot_on,
                    figtype=figtype,
                    methods=methods_ext,
                    methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                    phes=phes_ext,
                    phed=pconfig.phed(phes_ext, 'capital2'),
                    col_wrap=1,
                    palette=pconfig.palette_hls_2_methods(),
                    sharey=None,
                    ylabel=pconfig.accd()[plot_on],
                    xtick_label_format='norot',
                    xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                    ytick_ntick=4,
                    legend='horizontal',
                    legend_ncol=1)
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / (plot_on + '.integ-add-nonadd_phes3.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

    if 'poster' in pformats:
        # figtypes = ['2-4-1.2']
        # figtypes = ['3-4-1']
        figtypes = ['3-3.5-1.5']
        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                    'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                methods_ext = ['boosting_integrate', 'boosting_nonadd', 'boosting_add']
                if set(methods_ext) <= set(methods):
                    phes_ext = ['ra', 'psr', 'gout']
                    fig = plot_acc.plot_acc_box_strip_cat(
                        stats,
                        plot_on,
                        figtype=figtype,
                        methods=methods_ext,
                        methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                        phes=phes_ext,
                        phed=pconfig.phed(phes_ext, 'capital2'),
                        col_wrap=1,
                        palette=pconfig.palette_hls_2_methods(),
                        sharey=None,
                        ylabel=pconfig.accd()[plot_on],
                        xtick_label_format='norot',
                        xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                        ytick_ntick=3,
                        legend='horizontal',
                        legend_ncol=1)
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-add-nonadd_phes3.poster.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

    if 'slide' in pformats:
        # figtypes = ['3-3.5-1.5']
        # figtypes = ['4-3.5-1']
        figtypes = ['4-4-0.8']
        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                    'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                methods_ext = ['boosting_integrate', 'boosting_nonadd', 'boosting_add']
                if set(methods_ext) <= set(methods):
                    phes_ext = ['ra', 'psr', 'gout']
                    fig = plot_acc.plot_acc_box_strip_cat(
                        stats,
                        plot_on,
                        figtype=figtype,
                        methods=methods_ext,
                        methodd=pconfig.methodd(methods_ext, 'capital'),
                        phes=phes_ext,
                        phed=pconfig.phed(phes_ext, 'capital3'),
                        # col_wrap=1,
                        col_wrap=len(phes_ext),
                        palette=pconfig.palette_hls_2_methods(),
                        sharey=None,
                        ylabel=pconfig.accd()[plot_on],
                        xtick_format='none',
                        xtick_label_format='none',
                        # xtick_label_format='norot',
                        # xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                        ytick_ntick=3,
                        legend='horizontal')
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-add-nonadd_phes3.slide.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

    # for chr6
    # integ, prev; 3 phes
    if 'poster' in pformats:
        if ('ex-chr6' in list(methodd_type.values())) or ('noneur' in list(methodd_type.values())):
            # figtypes = ['2-4-1.2']
            # figtypes = ['3-4-1.4']
            # figtypes = ['4-4-1.5']
            figtypes = ['5-3.5-1.5']
            for figtype in figtypes:
                with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                        'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                    methods_ext = methods
                    # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_integrate_ex-chr6'])
                    phes_ext = ['ra', 'psr', 'gout']
                    fig = plot_acc.plot_acc_box_strip_cat(
                        stats, plot_on, figtype=figtype,
                        methods=methods_ext,
                        methodd=pconfig.methodd(methods_ext, 'capital'),
                        phes=phes_ext,
                        phed=pconfig.phed(phes_ext, 'capital2'),
                        col_wrap=1,
                        palette=pconfig.palette_hls_2_methods(),
                        sharey=None,
                        ylabel=pconfig.accd()[plot_on],
                        xtick_format='none',
                        xtick_label_format='none',
                        # xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                        ytick_ntick=3,
                        legend='horizontal',
                        legend_ncol=1,
                        vlines=[len(methods) / 2 - 0.5],
                    )
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-prev_phes3.poster.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

    # for chr6
    # integ, prev; 3 phes
    if 'slide' in pformats:
        if ('ex-chr6' in list(methodd_type.values())) or ('noneur' in list(methodd_type.values())):
            # figtypes = ['2-4-1.2']
            # figtypes = ['3-4-1.4']
            # figtypes = ['4-4-1.5']
            # figtypes = ['5-3.5-1.5']
            figtypes = ['6-3.5-1.3']
            for figtype in figtypes:
                with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                        'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                    methods_ext = methods
                    # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_integrate_ex-chr6'])
                    phes_ext = ['ra', 'psr', 'gout', 'ibd']
                    fig = plot_acc.plot_acc_box_strip_cat(
                        stats, plot_on, figtype=figtype,
                        methods=methods_ext,
                        methodd=pconfig.methodd(methods_ext, 'capital'),
                        phes=phes_ext,
                        phed=pconfig.phed(phes_ext, 'capital2'),
                        col_wrap=2,
                        # col_wrap=len(phes_ext),
                        # col_wrap=1,
                        palette=pconfig.palette_hls_2_methods(),
                        sharey=None,
                        ylabel=pconfig.accd()[plot_on],
                        xtick_format='none',
                        xtick_label_format='none',
                        # xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                        ytick_ntick=3,
                        legend='horizontal',
                        # legend_ncol=1,
                        vlines=[len(methods) / 2 - 0.5],
                    )
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-prev_phes3.slide.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

            # for legend
            figtypes = ['5-3.5-1.5']
            for figtype in figtypes:
                with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'xtick.labelsize': 18, 'legend.fontsize': 24,
                                                        'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                    methods_ext = methods[:int(len(methods) / 2)]
                    # methods_ext = methods
                    # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_integrate_ex-chr6'])
                    phes_ext = ['ra', 'psr', 'gout']
                    fig = plot_acc.plot_acc_box_strip_cat(
                        stats, plot_on, figtype=figtype,
                        methods=methods_ext,
                        methodd=pconfig.methodd(methods_ext, 'capital'),
                        phes=phes_ext,
                        phed=pconfig.phed(phes_ext, 'capital2'),
                        col_wrap=len(phes_ext),
                        # col_wrap=1,
                        palette=pconfig.palette_hls_2_methods(),
                        sharey=None,
                        ylabel=pconfig.accd()[plot_on],
                        xtick_format='none',
                        xtick_label_format='none',
                        # xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                        ytick_ntick=3,
                        legend='horizontal',
                        # legend_ncol=1,
                        # vlines=[len(methods) / 2 - 0.5],
                    )
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.integ-prev_phes3_legend.slide.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()


# def validate_acc_simu(droot, pformats, stats, plot_on, methods, phes):
def validate_acc_simu(droot, pformats, stats, plot_on, methods):

    phes = list(stats['phe'].unique())

    # all
    if 'paper' in pformats:
        # width=180mm=7inch
        # 5inch(default) * 6 =30inch
        # scale=7/30=1/4.3
        #
        # minimum font ytick.labelsize=18
        # 18pt*1/4.3=4.2pt
        # shold make this size >5pt as instructed by Nature
        # -> but looks fine in pptx
        #

        # FIXME: get legend of methods here

        figtype = '3-5-0.8'
        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            # stats_h2 = stats.loc[stats['h2'] == h2, :]
            fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
                                                  methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
                                                  methodd=pconfig.methodd(methods, 'capital_roman'),
                                                  phes=phes[:3],
                                                  phed=dict(zip(phes, phes)),
                                                  # phed=pconfig.phed(phes, 'capital2'),
                                                  # col_wrap=6,
                                                  palette=pconfig.palette_hls_2_methods(),
                                                  sharey=None,
                                                  ylabel=pconfig.accd()[plot_on],
                                                  xtick_label_format='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                  xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
                                                  ytick_ntick=4,
                                                  # ylim=(0.0, 0.25),
                                                  # legend=True,
                                                  legend='horizontal',
                                                  )
            # remove marker
            for item in fig._legend.legendHandles:
                item.set_visible(False)
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / (plot_on + '.legend.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

        # if len(methods) >= 10:
        #    figtypes = ['3-5-0.8', '4-5-1.5']
        # else:

        figtypes = ['3-5-1']
        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                    'xtick.major.width': 1.5, 'xtick.major.size': 8}):
                h2s = list(stats['h2'].unique())
                for h2 in h2s:
                    stats_h2 = stats.loc[stats['h2'] == h2, :]
                    fig = plot_acc.plot_acc_simu_box_strip_cat(stats_h2, plot_on, figtype=figtype,
                                                               methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
                                                               methodd=pconfig.methodd(methods, 'capital2_roman'),
                                                               phes=phes,
                                                               phed=dict(zip(phes, phes)),
                                                               # phed=pconfig.phed(phes, 'capital2'),
                                                               # col_wrap=6,
                                                               palette=pconfig.palette_hls_2_methods(),
                                                               sharey=None,
                                                               ylabel=pconfig.accd()[plot_on],
                                                               xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                               xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
                                                               ytick_ntick=4,
                                                               # ylim=(0.0, 0.25),
                                                               ylim=(0.0, float(h2) * 1.05),
                                                               #hlines=[float(h2)], # since no rationale for nagelkerke
                                                               # legend=True,
                                                               legend='horizontal',
                                                               )
                    # 0.8 time fontheight
                    fig.tight_layout(h_pad=0.8, w_pad=0)

                    fname = droot / (plot_on + '.all.h2-' + str(h2) + '.paper.figtype' + figtype)
                    logger.debug('fplot: {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()

            # col_wrap
            # with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
            #                                        'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            #    fig = plot_acc.plot_acc_box_strip_cat(stats, plot_on, figtype=figtype,
            #                                          methods=methods,  # methodd=pconfig.methodd(methods, 'capital'),
            #                                          methodd=pconfig.methodd(methods, 'capital2_roman'),
            #                                          phes=phes, phed=pconfig.phed(phes, 'capital2'),
            #                                          col_wrap=min(4, len(phes)),
            #                                          palette=pconfig.palette_hls_2_methods(),
            #                                          sharey=None,
            #                                          ylabel=pconfig.accd()[plot_on],
            #                                          xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
            #                                          xtick_label=pconfig.create_xtick_labels(methods, 'roman_number'),
            #                                          ytick_ntick=4,
            #                                          legend='horizontal',
            #                                          legend_ncol=math.ceil(len(methods) / 2))
            #    # 0.8 time fontheight
            #    fig.tight_layout(h_pad=0.8, w_pad=0)

            #    fname = droot / (plot_on + '.all.paper.col_wrap-4.figtype' + figtype)
            #    logger.debug('fplot: {}'.format(fname))
            #    pio.savefig_png_svg(fig, fname)
            #    plt.close()


def create_hlines_odds_ratio():
    hmax, hmin = 10.0, 0.0
    hlines = list(np.linspace(hmin, hmax, num=int((hmax - hmin) / 0.5) + 1))
    return hlines


# pformats=['paper', 'slide', 'poster']
def validate_odds_ratio(droot, pformats, stats, methods, phes):
    # TODO: avoid using margin_title_font_prop

    # all
    if 'paper' in pformats:
        figtype = '5-5-1.5'
        with sns.plotting_context(context='poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                        'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            hlines = create_hlines_odds_ratio()
            # validate_paper:3693
            fig = plot_odds_ratio.plot_odds_ratio_box_strip_cat(stats, figtype=figtype,
                                                                methods=methods, methodd=pconfig.methodd(methods, kind='capital'),
                                                                phes=phes, phed=pconfig.phed(phes, kind='capital2'),
                                                                props='1,3,5,10'.split(','),
                                                                col_wrap=4,
                                                                sharey=None,
                                                                palette=pconfig.palette_hls_2_methods(),
                                                                bottom_label='Genetic liability threshold in population stratification',
                                                                xlabel='',  # xlabel=' ',
                                                                ylabel='Odds ratio',
                                                                xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                # xtick_label=pconfig.create_xtick_labels(n=len(methods), kind='roman_number'),
                                                                xtick_label=pconfig.create_xtick_labels(methods, kind='roman_number'),
                                                                ytick_ntick=4,
                                                                legend='horizontal',
                                                                hlines=hlines,
                                                                margin_titile_font_prop=2,
                                                                )
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('or.all.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    # integ, add, nonadd+prev, 12phe
    if 'paper' in pformats:
        for figtype in ['9-6-2', '10-6-2.2']:
            # with sns.plotting_context(context='poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
            #                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            with sns.plotting_context(context='poster', rc={'font.size': 22,  # for suptitle used as xtick label
                                                            }):

                hlines = create_hlines_odds_ratio()
                methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
                fig = plot_odds_ratio.plot_odds_ratio_box_strip_cat(stats, figtype=figtype,
                                                                    methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital2_roman'),
                                                                    phes=pconfig.phes_12(), phed=pconfig.phed(pconfig.phes_12(), 'capital2'),
                                                                    props='1,3,5,10'.split(','),
                                                                    col_wrap=2,
                                                                    palette=pconfig.palette_hls_2_methods(),
                                                                    sharey=None,
                                                                    bottom_label='Genetic liability threshold in population stratification',
                                                                    xlabel='',  # xlabel=' ',
                                                                    ylabel='Odds ratio',
                                                                    xtick='norot',
                                                                    # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                    # xtick_label=pconfig.create_xtick_labels(n=len(methods), kind='roman_number'),
                                                                    xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                                    ytick_ntick=4,
                                                                    legend='horizontal',
                                                                    hlines=hlines,
                                                                    margin_titile_font_prop=2,
                                                                    )
                fig.tight_layout(rect=(0, 0.02, 0.95, 1))
                # 0.8 time fontheight
                # fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / ('or.integ-add-nonadd-prev_12phe.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

    # integ, add, nonadd+prev, phes3
    if 'paper' in pformats:
        figtype = '14-4.5-3'
        # with sns.plotting_context(context='poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
        #                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
        with sns.plotting_context(context='poster'):
            hlines = create_hlines_odds_ratio()
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'atm']
            fig = plot_odds_ratio.plot_odds_ratio_box_strip_cat(stats, figtype=figtype,
                                                                methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                                phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                                props='1,3,5,10'.split(','),
                                                                col_wrap=1,
                                                                palette=pconfig.palette_hls_2_methods(),
                                                                sharey=None,
                                                                xlabel='Genetic liability threshold in population stratification',
                                                                ylabel='Odds ratio',
                                                                xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                                ytick_ntick=4,
                                                                legend=False,
                                                                # legend='horizontal',
                                                                hlines=hlines,
                                                                margin_titile_font_prop=2,
                                                                )
            fig.tight_layout(rect=(0, 0.02, 0.95, 1))
            # 0.8 time fontheight
            # fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('or.integ-add-nonadd-prev_phes2.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    # integ, add, nonadd+prev, phes3
    if 'paper' in pformats:
        figtype = '14-4.5-3'
        # with sns.plotting_context(context='poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
        #                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
        with sns.plotting_context(context='poster'):
            hlines = create_hlines_odds_ratio()
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'gout']
            fig = plot_odds_ratio.plot_odds_ratio_box_strip_cat(stats, figtype=figtype,
                                                                methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                                phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                                props='1,3,5,10'.split(','),
                                                                col_wrap=1,
                                                                palette=pconfig.palette_hls_2_methods(),
                                                                sharey=None,
                                                                xlabel='Genetic liability threshold in population stratification',
                                                                ylabel='Odds ratio',
                                                                xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                                ytick_ntick=4,
                                                                legend=False,
                                                                # legend='horizontal',
                                                                hlines=hlines,
                                                                margin_titile_font_prop=2,
                                                                )
            fig.tight_layout(rect=(0, 0.02, 0.95, 1))
            # 0.8 time fontheight
            # fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('or.integ-add-nonadd-prev_phes3.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    # integ + prev, phes2
    if 'paper' in pformats:
        figtype = '14-4.5-3'
        # with sns.plotting_context(context='poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
        #                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
        with sns.plotting_context(context='poster'):
            hlines = create_hlines_odds_ratio()
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'atm']
            fig = plot_odds_ratio.plot_odds_ratio_box_strip_cat(stats, figtype=figtype,
                                                                methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                                phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                                props='1,3,5,10'.split(','),
                                                                col_wrap=1,
                                                                palette=pconfig.palette_hls_2_methods(),
                                                                sharey=None,
                                                                xlabel='Genetic liability threshold in population stratification',
                                                                ylabel='Odds ratio',
                                                                xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                                ytick_ntick=4,
                                                                legend=False,
                                                                # legend='horizontal',
                                                                hlines=hlines,
                                                                margin_titile_font_prop=2,
                                                                )
            fig.tight_layout(rect=(0, 0.02, 0.95, 1))
            # 0.8 time fontheight
            # fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('or.integ-prev_phes2.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    # integ + prev, phes3
    if 'paper' in pformats:
        figtype = '14-4.5-3'
        # with sns.plotting_context(context='poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
        #                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
        with sns.plotting_context(context='poster'):
            hlines = create_hlines_odds_ratio()
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'gout']
            fig = plot_odds_ratio.plot_odds_ratio_box_strip_cat(stats, figtype=figtype,
                                                                methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                                phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                                props='1,3,5,10'.split(','),
                                                                col_wrap=1,
                                                                palette=pconfig.palette_hls_2_methods(),
                                                                sharey=None,
                                                                xlabel='Genetic liability threshold in population stratification',
                                                                ylabel='Odds ratio',
                                                                xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                                xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                                ytick_ntick=4,
                                                                legend=False,
                                                                # legend='horizontal',
                                                                hlines=hlines,
                                                                margin_titile_font_prop=2,
                                                                )
            fig.tight_layout(rect=(0, 0.02, 0.95, 1))
            # 0.8 time fontheight
            # fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('or.integ-prev_phes3.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()


def create_hlines_nsnv():
    hmax = 9
    hmin = 1
    hline_ys = list(np.logspace(hmin, hmax, num=(hmax - hmin + 1), endpoint=True))
    return hline_ys


def validate_nsnv(droot, pformats, stats, methods, phes):

    # all
    if 'paper' in pformats:
        figtype = '5-5-1.0'
        with sns.plotting_context('poster', rc={'ytick.labelsize': 18, 'legend.fontsize': 24,
                                                'xtick.major.width': 1.5, 'xtick.major.size': 8}):
            fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
                                                    methods=methods, methodd=pconfig.methodd(methods, 'capital2_roman'),
                                                    phes=phes, phed=pconfig.phed(phes, 'capital2'),
                                                    col_wrap=min(6, len(phes)),
                                                    # col_wrap=6,
                                                    palette=pconfig.palette_hls_2_methods(),
                                                    sharey=True,
                                                    ylabel='# SNVs',
                                                    xtick_label_format='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number(len(methods)),
                                                    xtick_label=pconfig.create_xtick_labels(methods, kind='roman_number'),
                                                    legend='horizontal')
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('nsnv.all.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    if 'paper' in pformats:
        figtypes = ['5-5-1.0', '4-5-1.5']
        for figtype in figtypes:
            with sns.plotting_context('poster', rc={'xtick.labelsize': 18}):
                # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
                fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
                                                        methods=methods, methodd=pconfig.methodd(methods, kind='capital2_roman'),
                                                        phes=phes, phed=pconfig.phed(phes, kind='capital2'),
                                                        col_wrap=min(6, len(phes)),
                                                        # col_wrap=6,
                                                        palette=pconfig.palette_hls_2_methods(),
                                                        palette_strip=pconfig.palette_hls_2_methods_strip_nsnv(),
                                                        sharey=True, ylog=True,
                                                        ylabel='# SNVs',
                                                        xtick_label_format='norot',
                                                        # xtick_label=pconfig.create_xtick_labels_roman_number_boosting(methods_ext),
                                                        xtick_label=pconfig.create_xtick_labels(methods=methods, kind='roman_number'),
                                                        ytick=[100, 10000, 1000000],
                                                        ytick_label=['$10^2$', '$10^4$', '$10^6$'],
                                                        legend='horizontal',
                                                        hlines=create_hlines_nsnv(),
                                                        vlines=[2.5])
                # 0.8 time fontheight
                fig.tight_layout(h_pad=0.8, w_pad=0)

                fname = droot / ('nsnv.all.log.paper.figtype' + figtype)
                logger.debug('fplot: {}'.format(fname))
                pio.savefig_png_svg(fig, fname)
                plt.close()

    # logger.debug('tmp integ_add_nonadd_prev')
    # integ, add, nonadd + previous methods, 12phe
    # if 'paper' in pformats:
    #    figtype = '5-5-1.0'
    #    with sns.plotting_context('poster', rc={'xtick.labelsize': 18}):
    #        methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
    #        fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
    #                                                methods=methods_ext, methodd=pconfig.methodd(methods_ext, kind='capital_roman'),
    #                                                phes=pconfig.phes_12(), phed=pconfig.phed(pconfig.phes_12(), kind='capital2'),
    #                                                col_wrap=6,
    #                                                palette=pconfig.palette_hls_2_methods(),
    #                                                palette_strip=pconfig.palette_hls_2_methods_strip_nsnv(),
    #                                                sharey=True, ylog=True,
    #                                                ylabel='# SNVs',
    #                                                xtick='norot',  # xtick_label=pconfig.create_xtick_labels_roman_number_boosting(methods_ext),
    #                                                xtick_label=pconfig.create_xtick_labels(methods=methods_ext, kind='roman_number'),
    #                                                ytick=[100, 10000, 1000000],
    #                                                ytick_label=['$10^2$', '$10^4$', '$10^6$'],
    #                                                legend='horizontal',
    #                                                hlines=create_hlines_nsnv(),
    #                                                vlines=[2.5])
    #        # 0.8 time fontheight
    #        fig.tight_layout(h_pad=0.8, w_pad=0)

    #        fname = droot / ('nsnv.integ-add-nonadd-prev_12phe.paper.figtype' + figtype)
    #        logger.debug('fplot: {}'.format(fname))
    #        pio.savefig_png_svg(fig, fname)
    #        plt.close()

    # integ + prev, phes2
    if 'paper' in pformats:
        figtype = '10-4.5-1.3'
        with sns.plotting_context('poster', rc={'xtick.labelsize': 18}):
            # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'atm']
            fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
                                                    methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                    phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                    col_wrap=1,
                                                    palette=pconfig.palette_hls_2_methods(),
                                                    palette_strip=pconfig.palette_hls_2_methods_strip_nsnv(),
                                                    sharey=True, ylog=True,
                                                    ylabel='# SNVs',
                                                    xtick_label_format='norot',
                                                    # xtick_label=pconfig.create_xtick_labels_roman_number_boosting(methods_ext),
                                                    xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                    ytick=[100, 10000, 1000000],
                                                    ytick_label=['$10^2$', '$10^4$', '$10^6$'],
                                                    legend=False,
                                                    # legend='horizontal',
                                                    hlines=create_hlines_nsnv(),
                                                    # vlines=[2.5]
                                                    )
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('nsnv.integ-prev_phes1.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    # integ + prev, phes2
    if 'paper' in pformats:
        figtype = '10-4.5-1.3'
        with sns.plotting_context('poster', rc={'xtick.labelsize': 18}):
            # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'gout']
            fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
                                                    methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                    phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                    col_wrap=1,
                                                    palette=pconfig.palette_hls_2_methods(),
                                                    palette_strip=pconfig.palette_hls_2_methods_strip_nsnv(),
                                                    sharey=True, ylog=True,
                                                    ylabel='# SNVs',
                                                    xtick_label_format='norot',
                                                    # xtick_label=pconfig.create_xtick_labels_roman_number_boosting(methods_ext),
                                                    xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                    ytick=[100, 10000, 1000000],
                                                    ytick_label=['$10^2$', '$10^4$', '$10^6$'],
                                                    legend=False,
                                                    # legend='horizontal',
                                                    hlines=create_hlines_nsnv(),
                                                    # vlines=[2.5]
                                                    )
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('nsnv.integ-prev_phes3.paper.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    if 'poster' in pformats:
        # figtype = '1-4.5-1.3'
        # figtype = '2-4-1.5'
        # figtype = '3-4-1.2'
        # figtype = '4-3.5-1.3'
        figtype = '4-3.5-1.5'
        with sns.plotting_context('poster', rc={'xtick.labelsize': 18}):
            # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'gout']
            fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
                                                    methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital_roman'),
                                                    phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                    col_wrap=1,
                                                    palette=pconfig.palette_hls_2_methods(),
                                                    palette_strip=pconfig.palette_hls_2_methods_strip_nsnv(),
                                                    sharey=True, ylog=True,
                                                    ylabel='# SNVs',
                                                    xtick_format='none',
                                                    xtick_label_format='none',
                                                    # xtick_label=pconfig.create_xtick_labels_roman_number_boosting(methods_ext),
                                                    # xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                    ytick=[100, 10000, 1000000],
                                                    ytick_label=['$10^2$', '$10^4$', '$10^6$'],
                                                    legend=False,
                                                    # legend='horizontal',
                                                    hlines=create_hlines_nsnv(),
                                                    # vlines=[2.5]
                                                    )
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('nsnv.integ-prev_phes3.poster.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()

    if 'slide' in pformats:
        # figtype = '1-4.5-1.3'
        # figtype = '2-4-1.5'
        # figtype = '3-4-1.2'
        # figtype = '4-3.5-1.3'
        figtype = '4-3.5-1.5'
        with sns.plotting_context('poster', rc={'xtick.labelsize': 18}):
            # methods_ext = extract_boosing_methods(methods, ['boosting_integrate', 'boosting_add', 'boosting_nonadd'])
            methods_ext = extract_boosing_methods(methods, ['boosting_integrate'])
            # phes_ext = ['ra', 'atm', 'psr']
            phes_ext = ['ra', 'psr', 'gout', 'atm']
            fig = plot_nsnv.plot_nsnv_box_strip_cat(stats, figtype=figtype,
                                                    methods=methods_ext, methodd=pconfig.methodd(methods_ext, 'capital'),
                                                    phes=phes_ext, phed=pconfig.phed(phes_ext, 'capital2'),
                                                    col_wrap=2,
                                                    palette=pconfig.palette_hls_2_methods(),
                                                    palette_strip=pconfig.palette_hls_2_methods_strip_nsnv(),
                                                    sharey=True, ylog=True,
                                                    ylabel='# SNVs',
                                                    xtick_format='none',
                                                    xtick_label_format='none',
                                                    # xtick_label=pconfig.create_xtick_labels_roman_number_boosting(methods_ext),
                                                    # xtick_label=pconfig.create_xtick_labels(methods_ext, 'roman_number'),
                                                    ytick=[100, 10000, 1000000],
                                                    ytick_label=['$10^2$', '$10^4$', '$10^6$'],
                                                    legend=False,
                                                    # legend='horizontal',
                                                    hlines=create_hlines_nsnv(),
                                                    # vlines=[2.5]
                                                    )
            # 0.8 time fontheight
            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('nsnv.integ-prev_phes3.slide.figtype' + figtype)
            logger.debug('fplot: {}'.format(fname))
            pio.savefig_png_svg(fig, fname)
            plt.close()


def validate_genetic_model(droot, pformats, stats, phes):

    # TMP cvi=0 only
    gmodels = ['overrec', 'rec', 'add', 'dom', 'overdom']

    def stats_gmodels(stats, cvi):
        logger.debug('stats {}'.format(stats))
        logger.debug('stats {}'.format(stats.info()))
        stats_cvi = stats.loc[stats['cvi'] == cvi, ['phe'] + gmodels]
        # check if only one stat for each phe
        assert (~stats_cvi['phe'].duplicated()).all()
        gmodel = stats_cvi.set_index('phe')
        return gmodel

    cvi = 0
    genetic_models = stats_gmodels(stats, cvi)
    # genetic_models = stats

    phes_ext = pconfig.phes_12()
    genetic_models = genetic_models.loc[phes_ext, :]

    # all
    if 'paper' in pformats:
        # """
        # hatch on barplot
        # https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
        # https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
        # """

        figtype = '1-8-1'
        fig = plot_genetic_model.genetic_model_prop(genetic_models, figtype)
        fig.tight_layout()

        fname = droot / ('genetic-model.12phe.percent_bar.nolegend.paper.gwasmethods.figtype' + figtype)
        pio.savefig_png_svg(fig, fname)
        plt.close()

        # TODO: legend

        # ['overdom', 'dom','add','rec', 'overrec']
        # should be different colors as odds ratio
        # mypalette4 = [
        #    cb[6], cb[4], cb[2], cb[0], cb[9]
        # ]
        # for stacked bar
        # mypalette4r = list(reversed(mypalette4))

        # first, bottom
        # change legend order later
        # columns_order = ['overrec', 'rec', 'add', 'dom', 'overdom']
        # columns_order = ['overdom', 'dom', 'add', 'rec', 'overrec']
        # columns_order = ['add', 'overdom', 'dom', 'overrec', 'rec']
        # genetic_models = genetic_models.loc[:, columns_order]
        # genetic_models.index = genetic_models.index.map(pconfig._phed_capital3())
        # genetic_models.index = genetic_models.index.map(get_phed())
        # genetic_models.columns = genetic_models.columns.map(pconfig.genetic_modeld())
        # hatchs = ['//', '\\\\', '', '||', '--', ]
        # hatchs = ['//', '\\\\', '..', '||', '--', ]
        # hatchs = ['--', '++', '.', '//', '\\']
        # order is inversed
        # mypalette3_alpha = [v + (1.0,) for v in mypalette4r]
        # color_hatch = dict(zip(mypalette3_alpha, hatchs))
        # with sns.plotting_context('poster'):
        #    with sns.color_palette(mypalette4r, 5):
        #        print('gm', genetic_models)
        #        genetic_models_percent = genetic_models.copy()
        #        for index in genetic_models_percent.index:
        #            print('index', index)
        #            print('row', genetic_models_percent.loc[index, :])
        #            genetic_models_percent.loc[index, :] = genetic_models_percent.loc[index, :] * 100 / genetic_models_percent.loc[index, :].sum()
        #        print('genetic_models_percent', genetic_models_percent)

        #        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        #        # edgecolor is for hatch color but also change edge
        #        ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False, edgecolor='w')
        #        # ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False)
        #        for mi, bar in enumerate(ax.patches):
        #            hatch = color_hatch[bar.get_facecolor()]
        #            bar.set_hatch(hatch)
        #        ax.set(xlabel='Phenotype',
        #               ylabel='Proportion of SNVs')

        #        ax.set_ylim((0, 100))
        #        ax.set_yticks((0, 25, 50, 75, 100))
        #        ax.set_yticklabels(('0%', '25%', '50%', '75%', '100%'))
        #        # https://analytics-note.xyz/programming/matplotlib-tick-rotation/
        #        ax.yaxis.set_tick_params(labelsize=20)
        #        ax.xaxis.set_tick_params(labelsize=15, rotation=90)
        #        fig.tight_layout()
        #        # fig.savefig(fout + '.percent.bar.nolegend.xrot.png')
        #        # fig.savefig(fout + '.percent.bar.nolegend.xrot.pdf')

        #        figtype = '1-8-1'

        #        fname = droot / ('genetic-model.12phe.percent_bar.nolegend.paper.gwasmethods.figtype' + figtype)
        #        pio.savefig_png_svg(fig, fname)
        #        plt.close()

    if 'poster' in pformats:

        figtype = '1-8-1'
        fig = plot_genetic_model.genetic_model_prop(genetic_models, figtype, xlabel_fontsize=28, ylabel_fontsize=28)
        fig.tight_layout()

        fname = droot / ('genetic-model.12phe.percent_bar.nolegend.poster.gwasmethods.figtype' + figtype)
        pio.savefig_png_svg(fig, fname)
        plt.close()

        # """
        # hatch on barplot
        # https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
        # https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
        # """

        # cb = sns.color_palette("colorblind", 10)
        # ['overdom', 'dom','add','rec', 'overrec']
        # should be different colors as odds ratio
        # mypalette4 = [
        #    cb[6], cb[4], cb[2], cb[0], cb[9]
        # ]
        # for stacked bar
        # mypalette4r = list(reversed(mypalette4))

        # first, bottom
        # change legend order later
        # columns_order = ['overrec', 'rec', 'add', 'dom', 'overdom']
        # columns_order = ['overdom', 'dom', 'add', 'rec', 'overrec']
        # columns_order = ['add', 'overdom', 'dom', 'overrec', 'rec']
        # genetic_models = genetic_models.loc[:, columns_order]
        # genetic_models.index = genetic_models.index.map(pconfig._phed_capital3())
        # genetic_models.index = genetic_models.index.map(get_phed())
        # genetic_models.columns = genetic_models.columns.map(pconfig.genetic_modeld())
        # hatchs = ['//', '\\\\', '', '||', '--', ]
        # hatchs = ['//', '\\\\', '..', '||', '--', ]
        # hatchs = ['--', '++', '.', '//', '\\']
        # order is inversed
        # mypalette3_alpha = [v + (1.0,) for v in mypalette4r]
        # color_hatch = dict(zip(mypalette3_alpha, hatchs))
        # with sns.plotting_context('poster'):
        #    with sns.color_palette(mypalette4r, 5):
        #        print('gm', genetic_models)
        #        genetic_models_percent = genetic_models.copy()
        #        for index in genetic_models_percent.index:
        #            print('index', index)
        #            print('row', genetic_models_percent.loc[index, :])
        #            genetic_models_percent.loc[index, :] = genetic_models_percent.loc[index, :] * 100 / genetic_models_percent.loc[index, :].sum()
        #        print('genetic_models_percent', genetic_models_percent)

        #        # TMP
        #        # figtype = '1-8-1'
        #        # fig, ax = plt.subplots(1, 1, figsize=(8, 8))

        #        # figtype = '1-7-0.7'
        #        # fig, ax = plt.subplots(1, 1, figsize=(5 * 1.4, 5))
        #        # figtype = '2-5-1.4'
        #        figtype = '1-8-1'
        #        fig, ax = plt.subplots(1, 1, figsize=(8, 8))

        #        # edgecolor is for hatch color but also change edge
        #        ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False, edgecolor='w')
        #        # ax = genetic_models_percent.plot.bar(stacked=True, ax=ax, legend=False)
        #        for mi, bar in enumerate(ax.patches):
        #            hatch = color_hatch[bar.get_facecolor()]
        #            bar.set_hatch(hatch)
        #        # ax.set(xlabel='Phenotype',
        #            # ylabel='Proportion of SNVs')
        #        ax.set_xlabel('Phenotype', fontsize=28)
        #        ax.set_ylabel('Proportion of SNVs', fontsize=28)

        #        ax.set_ylim((0, 100))
        #        ax.set_yticks((0, 25, 50, 75, 100))
        #        ax.set_yticklabels(('0%', '25%', '50%', '75%', '100%'))
        #        # https://analytics-note.xyz/programming/matplotlib-tick-rotation/
        #        ax.yaxis.set_tick_params(labelsize=20)
        #        ax.xaxis.set_tick_params(labelsize=15, rotation=90)
        #        fig.tight_layout()
        #        # fig.savefig(fout + '.percent.bar.nolegend.xrot.png')
        #        # fig.savefig(fout + '.percent.bar.nolegend.xrot.pdf')

        #        fname = droot / ('genetic-model.12phe.percent_bar.nolegend.poster.gwasmethods.figtype' + figtype)
        #        pio.savefig_png_svg(fig, fname)
        #        plt.close()

    if 'slide' in pformats:

        figtype = '1-8-1'
        fig = plot_genetic_model.genetic_model_prop(genetic_models, figtype, xlabel_fontsize=28, ylabel_fontsize=28)
        fig.tight_layout()

        fname = droot / ('genetic-model.12phe.percent_bar.nolegend.slide.gwasmethods.figtype' + figtype)
        pio.savefig_png_svg(fig, fname)
        plt.close()


def validate_genetic_model_vs_nsnv(droot, pformats, stats, phes):

    # TMP cvi=0 only
    gmodels = ['overrec', 'rec', 'add', 'dom', 'overdom']

    def stats_gmodels(stats, cvi):
        logger.debug('stats {}'.format(stats))
        logger.debug('stats {}'.format(stats.info()))
        stats_cvi = stats.loc[stats['cvi'] == cvi, ['phe', 'cvi', 'nsnv_use'] + gmodels]
        # stats_cvi = stats.loc[stats['cvi'] == cvi, ['phe'] + gmodels]
        # check if only one stat for each phe
        assert (~stats_cvi['phe'].duplicated()).all()
        gmodel = stats_cvi.set_index('phe')
        return gmodel

    def nonadd_count(stats):
        stats.loc[:, 'nonadd'] = stats['overrec'] + stats['rec'] + stats['dom'] + stats['overdom']
        stats.loc[:, 'nonadd_prop'] = stats['nonadd'] / (stats['nonadd'] + stats['add'])
        return stats

    cvi = 0
    genetic_models = stats_gmodels(stats, cvi)
    genetic_models = nonadd_count(genetic_models)
    # genetic_models = stats

    phes_ext = pconfig.phes_12()
    genetic_models = genetic_models.loc[phes_ext, :]

    logger.debug("genetic_models: {}".format(genetic_models))

    # all
    if 'paper' in pformats:
        """
        # hatch on barplot
        https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
        https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
        """

        with sns.plotting_context('poster'):

            figtype = '1-8-1'
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))

            ax = sns.scatterplot(data=genetic_models, x='nsnv_use', y='nonadd_prop', ax=ax)

            # should set log before get_xlim
            ax.set(xscale="log")
            xlim = ax.get_xlim()
            logger.debug('xlim: {}'.format(xlim))

            # regression line for log
            x = np.log10(genetic_models['nsnv_use'].to_numpy())
            y = genetic_models['nonadd_prop'].to_numpy()
            logger.debug('x,y:{},{}'.format(x.shape, y.shape))
            reg = LinearRegression().fit(x.reshape((-1, 1)), y)
            coef = reg.coef_
            intercept = reg.intercept_
            line_x_log = [0, 10]
            line_y = [coef * x + intercept for x in line_x_log]
            line_x = [np.power(10, x) for x in line_x_log]
            cb = sns.color_palette("colorblind", 10)
            ax.plot(line_x, line_y, color=cb[7])

            ax.set_xlim(xlim)
            ax.set_ylim([0.0, 1.0])
            # ax.set(xscale="log")
            ax.set_xlabel('# SNVs')
            ax.set_ylabel('Proportion of non-additive SNVs')

            fig.tight_layout(h_pad=0.8, w_pad=0)

            fname = droot / ('relplot.12phe.paper.figtype' + figtype)
            pio.savefig_png_svg(fig, fname)
            plt.close()


def validate_manhattan(droot, pformats, ssd, methods, phe, gmodels=None):

    if gmodels is None:
        gmodels = list(ssd.keys())

    # all
    if 'paper' in pformats:
        # figtype = '1-10-15'
        figtype = '1-10-1.5'
        for gmodel in gmodels:
            ss = ssd[gmodel]
            xaxis_skip = [14, 16, 18, 20, 22]
            # no effect
            # with sns.plotting_context(font_scale=2):
            with sns.plotting_context('talk'):
                fig = plot_manhattan.plot_manhattan(ss, figtype=figtype,
                                                    ycol_method='effect_abs',
                                                    methods=methods, methodd=pconfig.methodd(methods, kind='capital2'),
                                                    title=pconfig.phed([phe], kind='capital2')[phe],
                                                    ylim='ymax_pos',
                                                    xaxis_skip=xaxis_skip,
                                                    )
                fig.tight_layout()

                fname = droot / ('manhattan.all.' + gmodel + '.paper.figtype' + figtype)
                logger.debug('fname {}'.format(fname))
                pio.savefig_png_svg(fig, fname, dpi=1200)
                plt.close()

        # close to gwasmethods
        for gmodel in gmodels:
            ss = ssd[gmodel]
            xaxis_skip = [14, 16, 18, 20, 22]
            with sns.plotting_context('talk'):
                fig = plot_manhattan.plot_manhattan(ss, figtype=figtype,
                                                    methods=methods, methodd=pconfig.methodd(methods, kind='capital2'),
                                                    ycol_method='effect',
                                                    xaxis_skip=xaxis_skip,
                                                    ylim='ymax_pos',
                                                    )
                fig.tight_layout()

                # fname = droot / ('manhattan.all_positive_effect.' + gmodel + '.paper.figtype' + figtype)
                fname = droot / ('manhattan.positive_effect_only.' + gmodel + '.paper.figtype' + figtype)
                logger.debug('fname {}'.format(fname))
                pio.savefig_png_svg(fig, fname, dpi=1200)
                plt.close()

        # close to gwasmethods
        for gmodel in gmodels:
            ss = ssd[gmodel]
            xaxis_skip = [14, 16, 18, 20, 22]
            with sns.plotting_context('talk'):
                fig = plot_manhattan.plot_manhattan(ss, figtype=figtype,
                                                    methods=methods, methodd=pconfig.methodd(methods, kind='capital2'),
                                                    ycol_method='effect',
                                                    xaxis_skip=xaxis_skip,
                                                    )
                fig.tight_layout()

                fname = droot / ('manhattan.all_effect.' + gmodel + '.paper.figtype' + figtype)
                logger.debug('fname {}'.format(fname))
                pio.savefig_png_svg(fig, fname, dpi=1200)
                plt.close()


def create_heatmap_table(ss_wgt, gwas_model_est_criteria, gmodels=None):

    if gmodels is None:
        gmodels = ['overdom', 'dom', 'add', 'rec', 'overrec']

    if gwas_model_est_criteria == 'gwas':
        gwas_model_est = 'model_est_gwas'
    elif gwas_model_est_criteria == 'aic':
        gwas_model_est = 'model_est_aic'

    data = ss_wgt[['boost_model_est', gwas_model_est]].copy()
    # data = ss_wgt[['boost_model_est', 'model_est_gwas']].copy()
    data.loc[:, 'dummy'] = 1
    # print('data', data.head())

    stat_heatmap = pd.pivot_table(data,
                                  values='dummy',
                                  index='boost_model_est', columns=gwas_model_est,
                                  # index='boost_model_est', columns='model_est_gwas',
                                  aggfunc='count')

    logger.debug('data pivot: {}'.format(data))

    def fill_col_index(data, gmodels):
        for gmodel in gmodels:
            # fill index
            if gmodel not in data.index:
                data.loc[gmodel, :] = np.nan
            # fill col
            if gmodel not in data:
                data.loc[:, gmodel] = np.nan
        return data
    stat_heatmap = fill_col_index(stat_heatmap, gmodels)

    stat_heatmap = stat_heatmap.loc[gmodels, gmodels]

    # FIXME: use categorical type for boost_model_est
    # then, these are not necessary
    stat_heatmap = stat_heatmap.fillna(0.0)
    stat_heatmap = stat_heatmap.astype(int)
    logger.debug('stat_heatmap pivot: {}'.format(stat_heatmap))
    return stat_heatmap


def validate_genetic_model_gwas(droot, pformats, stat, phe):
    """
    stat=wgt_ss
    # TMP cvi=0 only
    """

    ss_wgt = stat

    gwas_model_est_criteria = 'gwas'

    # gmodels = ['overdom', 'dom', 'add', 'rec', 'overrec']
    gmodels = pconfig.genetic_models()

    ss_wgt = ss_wgt.copy()
    stat_heatmap = create_heatmap_table(ss_wgt, gwas_model_est_criteria, gmodels)

    
    stat.to_csv(droot / ('stat.tsv'), na_rep='NaN', sep='\t', index=False)
    #stat_paper_public=public_stat.stat_to_paper_public(stat,['phe','cvi','nsnv_use','overdom','dom','add','rec','overrec'])
    #stat_paper_public.to_csv(droot / ( 'stat_public.tsv'), na_rep='NaN', sep='\t', index=False)


    # all
    if 'paper' in pformats:
        """
        # hatch on barplot
        https://www.tutorialspoint.com/how-do-i-plot-hatched-bars-using-pandas-and-matplotlib
        https://stackoverflow.com/questions/22833404/how-do-i-plot-hatched-bars-using-pandas
        """
        with sns.plotting_context('poster', font_scale=1.0,
                                  # rc=dict(font='Calibri')
                                  ):
            # with sns.set(font='Calibri'):
            # logger.debug('scatter heatmap')
            # gwas_model_est_criteria = 'gwas'
            # fig = plot_genetic_model.plot_genetic_model_gwas_heatmap_gwasmethods(stat, gwas_model_est_criteria)
            # fname = droot / ('heatmap' + '.paper.gwasmethods.figtype' + '0')
            # logger.debug('fname {}'.format(fname))
            # pio.savefig_png_svg(fig, fname, dpi=1200)
            # plt.close()

            gwas_model_est_criteria = 'gwas'
            fig = plot_genetic_model.plot_genetic_model_gwas_heatmap(stat_heatmap,
                                                                     gwas_model_est_criteria,
                                                                     # gmodels=gmodels,
                                                                     gmodeld=pconfig.genetic_modeld2(),
                                                                     palette=pconfig.palette_hls_2_genetic_model(),
                                                                     )
            fname = droot / ('heatmap' + '.paper.figtype' + '0')
            logger.debug('fname {}'.format(fname))
            pio.savefig_png_svg(fig, fname, dpi=1200)
            plt.close()


def get_figsize_locus_zoom(figtype, plot_gene, gm_n, ncols=1):
    # +1 for gene
    if plot_gene:
        # gm_n = len(genetic_models)
        if figtype == 1:
            figwidth = 8
            figsize = (figwidth * ncols, 5 * gm_n + 2.5)
            # [1.0] for blank padding axis
            height_ratios = [5] * gm_n + [1.0] + [2.5]
        elif figtype == 2:
            figwidth = 5
            figsize = (figwidth * ncols, 2 * gm_n + 1.0)
            height_ratios = [2] * gm_n + [1.0]
        elif figtype == 3:
            figwidth = 8
            # 1.0 for title
            figsize = (figwidth * ncols, 1.0 + 2.5 * gm_n + 2.5)
            height_ratios = [2.5] * gm_n + [2.5]
        elif figtype == 4:
            figwidth = 6
            # 1.0 for title
            figsize = (figwidth * ncols, 1.0 + 4 * gm_n + 2.5)
            height_ratios = [2.5] * gm_n + [2.5]
        else:
            raise RuntimeError()
        # +2 for gene and black padding axis bet. man and gene
        fig, axs = plt.subplots(nrows=gm_n + 2, ncols=ncols, figsize=figsize, sharex='col',
                                sharey='row', gridspec_kw={'height_ratios': height_ratios})
        # fig, axs = plt.subplots(nrows=gm_n + 1,ncols=ncols, figsize=figsize, sharex='col',sharey='row', gridspec_kw={'height_ratios': height_ratios})
        # fig, axs = plt.subplots(nrows=gm_n + 1,ncols=ncols, figsize=figsize,
        # sharex='col',sharey=False, gridspec_kw={'height_ratios': height_ratios})
    else:
        # gm_n = len(genetic_models)
        if figtype == 1:
            figwidth = 8
            figsize = (figwidth * ncols, 5 * gm_n)
        elif figtype == 2:
            figwidth = 5
            figsize = (figwidth * ncols, 2 * gm_n)
        elif figtype == 3:
            # 1.0 for title
            figwidth = 8
            figsize = (figwidth * ncols, 1.0 + 2.5 * gm_n)
        else:
            raise RuntimeError()

        fig, axs = plt.subplots(nrows=gm_n, ncols=ncols, figsize=figsize, sharex='col', sharey='row')

    return fig, axs, figwidth


def validate_locus_zoom(droot, pformats, ssd, methods, phe, local_ranges, gene_exon):

    # all
    if 'paper' in pformats:
        for rangei, local_range in enumerate(local_ranges):
            figtype = '1'
            with sns.plotting_context('poster'):
                if len(local_range) == 1:
                    local_range_ = local_range[0]
                    if len(local_range_) == 2:
                        chrompos, gmodels = local_range_
                        annot = None
                        plot_gene = True
                    elif len(local_range_) == 3:
                        chrompos, gmodels, annot = local_range_
                        plot_gene = True
                    elif len(local_range_) == 4:
                        chrompos, gmodels, annot, plot_gene = local_range_
                    else:
                        raise RuntimeError('wrong len')
                    
                    # TMP
                    #for figtype in ['1', '9']:
                    for figtype in ['9']:
                        fig = plot_manhattan.plot_locus_zoom_multi(
                            ssd,
                            chrompos,
                            figtype,
                            phe_title=pconfig.phed([phe], kind='capital2')[phe],
                            genetic_models=gmodels,
                            gene_exon=gene_exon,
                            plot_loss=True,
                            plot_gene=plot_gene,
                            # plot_gene=True,
                            annot=annot,
                        )

                        fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.bfr-layout.figtype' + figtype)
                        #fname = droot / ('locus-zoom.all.multi_cols.range-' + str(rangei) + '.paper.bfr-layout.figtype' + figtype)
                        logger.debug('fname {}'.format(fname))
                        pio.savefig_png_svg(fig, fname, dpi=1200)

                        fig.tight_layout(h_pad=0.0)
                        fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.figtype' + figtype)
                        logger.debug('fname {}'.format(fname))
                        pio.savefig_png_svg(fig, fname, dpi=1200)

                elif len(local_range) <= 3:
                    fig = plot_manhattan.plot_locus_zoom_multi_cols(
                        ssd,
                        phe_title=pconfig.phed([phe], kind='capital2')[phe],
                        chromposs=local_range,
                        gene_exon=gene_exon,
                        plot_loss=True,
                        plot_gene=True,
                        # annot=annot,
                    )

                    fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.bfr-layout.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname, dpi=1200)

                    # fig.subplots_adjust(hspace=0.05)
                    # fig.subplots_adjust(hspace=0.07)
                    # fig.subplots_adjust(hspace=0.15)
                    # this makes hspace bet. man too large
                    # fig.tight_layout()
                    fig.subplots_adjust(hspace=0.1)
                    # fig.subplots_adjust(hspace=0.1,wspace=0.3)
                    fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname, dpi=1200)

                elif len(local_range) <= 6:
                    logger.debug('(local_range): {}'.format(local_range))
                    logger.debug('len(local_range): {}'.format(len(local_range)))
                    fig, axs = plot_manhattan.plot_locus_zoom_multi_cols(
                        ssd,
                        phe_title=pconfig.phed([phe], kind='capital2')[phe],
                        chromposs=local_range,
                        gene_exon=gene_exon,
                        plot_loss=True,
                        plot_gene=True,
                        # annot=annot,
                    )

                    fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.bfr-layout.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname, dpi=1200)

                    # fig.subplots_adjust(hspace=0.05)
                    # fig.subplots_adjust(hspace=0.07)
                    # fig.subplots_adjust(hspace=0.15)
                    # this makes hspace bet. man too large
                    # fig.tight_layout()
                    fig.subplots_adjust(hspace=0.1)
                    # fig.subplots_adjust(hspace=0.1,wspace=0.3)
                    fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname, dpi=1200)

                    # split rows
                    # this retuns 1dim list -> how to return two dim?
                    # axs=fig.axes
                    # logger.debug('axs {}'.format(axs))
                    # logger.debug('axs {}'.format(axs.shape))
                    ylim_0 = axs[0, 0].get_ylim()[1]
                    ylim_1 = axs[1, 0].get_ylim()[1]
                    ylims = (ylim_0, ylim_1)

                    fig, _ = plot_manhattan.plot_locus_zoom_multi_cols(
                        ssd,
                        phe_title=pconfig.phed([phe], kind='capital2')[phe],
                        chromposs=local_range[:3],
                        gene_exon=gene_exon,
                        plot_loss=True,
                        plot_gene=True,
                        ylims=ylims
                    )

                    fig.subplots_adjust(hspace=0.1)
                    fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.sub0' + '.paper.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname, dpi=1200)

                    fig, _ = plot_manhattan.plot_locus_zoom_multi_cols(
                        ssd,
                        phe_title=pconfig.phed([phe], kind='capital2')[phe],
                        chromposs=local_range[3:],
                        gene_exon=gene_exon,
                        plot_loss=True,
                        plot_gene=True,
                        ylims=ylims
                    )

                    fig.subplots_adjust(hspace=0.1)
                    fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.sub1' + '.paper.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname, dpi=1200)
                else:
                    raise NotImplementedError('Too much local_range: ', len(local_range))

    if 'poster' in pformats:
        for rangei, local_range in enumerate(local_ranges):
            # figtypes = ['5', '6']
            figtypes = ['5']
            with sns.plotting_context('poster'):
                for figtype in figtypes:
                    if len(local_range) == 1:
                        local_range_ = local_range[0]
                        if len(local_range_) == 2:
                            chrompos, gmodels = local_range_
                            annot = None
                            plot_gene = True
                        elif len(local_range_) == 3:
                            chrompos, gmodels, annot = local_range_
                            plot_gene = True
                        elif len(local_range_) == 4:
                            chrompos, gmodels, annot, plot_gene = local_range_
                        else:
                            raise RuntimeError('wrong len')
                        fig = plot_manhattan.plot_locus_zoom_multi(
                            ssd,
                            chrompos,
                            figtype,
                            phe_title=pconfig.phed([phe], kind='capital2')[phe],
                            genetic_models=gmodels,
                            gene_exon=gene_exon,
                            plot_loss=True,
                            plot_gene=plot_gene,
                            # plot_gene=True,
                            annot=annot,
                        )

                        fname = droot / ('locus-zoom.all.multi_cols.range-' + str(rangei) + '.paper.bfr-layout.figtype' + figtype)
                        logger.debug('fname {}'.format(fname))
                        pio.savefig_png_svg(fig, fname, dpi=1200)

                        fig.tight_layout(h_pad=0.0)
                        fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.paper.figtype' + figtype)
                        logger.debug('fname {}'.format(fname))
                        pio.savefig_png_svg(fig, fname, dpi=1200)
                    else:
                        pass
                        # raise NotImplementedError('Too much local_range: ', len(local_range))

    if 'slide' in pformats:
        for rangei, local_range in enumerate(local_ranges):
            # figtype = '1'
            # figtype = '4'
            # for figtype in ['4', '5']:
            for figtype in ['6']:
                # no effect
                # with sns.plotting_context('poster', rc={'axes.labelsize': 60}):
                with sns.plotting_context('poster'):
                    if len(local_range) == 1:
                        local_range_ = local_range[0]
                        if len(local_range_) == 2:
                            chrompos, gmodels = local_range_
                            annot = None
                            plot_gene = True
                        elif len(local_range_) == 3:
                            chrompos, gmodels, annot = local_range_
                            plot_gene = True
                        elif len(local_range_) == 4:
                            chrompos, gmodels, annot, plot_gene = local_range_
                        else:
                            raise RuntimeError('wrong len')

                        fig = plot_manhattan.plot_locus_zoom_cols(
                            ssd,
                            chrompos,
                            figtype,
                            phe_title=pconfig.phed([phe], kind='capital2')[phe],
                            genetic_models=gmodels,
                            gene_exon=gene_exon,
                            plot_loss=True,
                            plot_gene=plot_gene,
                            # plot_gene=True,
                            annot=annot,
                            annot_pos='middle'
                        )

                        # fig = plot_manhattan.plot_locus_zoom_multi(
                        #    ssd,
                        #    chrompos,
                        #    figtype,
                        #    phe_title=pconfig.phed([phe], kind='capital2')[phe],
                        #    genetic_models=gmodels,
                        #    gene_exon=gene_exon,
                        #    plot_loss=True,
                        #    plot_gene=plot_gene,
                        #    # plot_gene=True,
                        #    annot=annot,
                        # )

                        fname = droot / ('locus-zoom.all.multi_cols.range-' + str(rangei) + '.slide.bfr-layout.figtype' + figtype)
                        logger.debug('fname {}'.format(fname))
                        pio.savefig_png_svg(fig, fname, dpi=1200)

                        fig.tight_layout(h_pad=0.0)
                        fname = droot / ('locus-zoom.all.range-' + str(rangei) + '.slide.figtype' + figtype)
                        logger.debug('fname {}'.format(fname))
                        pio.savefig_png_svg(fig, fname, dpi=1200)
                    else:
                        pass


def validate_phes_correlation(droot, pformats, stat, sex, stat_name):

    if 'paper' in pformats:
        with sns.plotting_context('poster',  # font_scale=1.0,
                                  rc={'font.size': 10,
                                      'ytick.labelsize': 12, 'xtick.labelsize': 12
                                      # 'ytick.labelsize': 18, 'xtick.labelsize': 18
                                      }):
            figtype = '0-10-1'
            # figtype = '0-12-1'
            fig = plot_phes_correlation.plot_stat_heatmap(stat,
                                                          figtype=figtype,
                                                          phed=pconfig.phed(stat.index.to_list(), 'capital3'),
                                                          # palette=pconfig.palette_hls_2_genetic_model(),
                                                          sex=sex,
                                                          stat_name=stat_name,
                                                          statd=pconfig.statd(),
                                                          )
            fig.tight_layout()
            # fig.tight_layout(h_pad=-0.8, w_pad=0)
            fname = droot / (stat_name + '.' + sex + '.' + 'heatmap.paper.figtype' + figtype)
            logger.debug('fname {}'.format(fname))
            pio.savefig_png_svg(fig, fname, dpi=1200)
            plt.close()


def validate_h2(droot, pformats, stat):

    if 'slide' in pformats:
        with sns.plotting_context('poster',  # font_scale=1.0,
                                  rc={
                                      # 'font.size': 10,
                                      # 'ytick.labelsize': 12, 'xtick.labelsize': 12
                                      # 'ytick.labelsize': 18, 'xtick.labelsize': 18
                                  }):
            figtype = '1-6-1'
            # figtype = '0-10-1'
            # figtype = '0-12-1'
            fig = plot_h2.plot_scatter(stat,
                                       figtype=figtype,
                                       phed=pconfig.phed(stat['phe'].to_list(), 'capital3'),
                                       )
            fig.tight_layout()
            # fig.tight_layout(h_pad=-0.8, w_pad=0)
            fname = droot / ('h2.scatter.slide.figtype' + figtype)
            logger.debug('fname {}'.format(fname))
            pio.savefig_png_svg(fig, fname, dpi=1200)
            plt.close()


def validate_sample_score(droot, pformats, scored, methods, phe):
    # create rs

    # TODO: for
    dataset = 'ts'

    rs = None
    for mi, method in enumerate(scored):
        score = scored[method][dataset]
        if mi == 0:
            rs = score
        else:
            assert len(rs) == len(score)
            score_tmp = score[['id', 'score']].rename(columns={'score': method + '_score'}).copy()
            rs = pd.merge(rs, score_tmp[['id', method + '_score']], on='id', how='outer')

    rs.loc[:, 'status'] = rs['status'].map({0: 'control', 1: 'case'})

    logger.debug('rs {}'.format(rs))

    # (x_method, y_method) = ('boosting_add', 'boosting_nonadd')
    # (x_method, y_method) = ( 'boosting_nonadd','boosting_add')
    xys = [
        ('boosting_nonadd', 'boosting_add'),
        ('boosting_add', 'boosting_nonadd'),
        ('boosting_nonadd', 'snpnet'),
    ]

    if 'paper' in pformats:
        for (x_method, y_method) in xys:
            for stat in ['count', 'proportion']:
                # figtype = '3-5-1.0'
                figtype = '3'
                with sns.plotting_context('poster',
                                          # rc={'lines.linewidth': 2.0}
                                          # rc={'ytick.labelsize': 18, 'legend.fontsize':24}
                                          ):
                    fig = plot_sample_score.plot_sample_score_summary(
                        rs,
                        x_method + '_score', y_method + '_score',
                        xlabel=pconfig.methodd(methods, kind='capital')[x_method] + '\nsample score',
                        ylabel=pconfig.methodd(methods, kind='capital')[y_method] + '\nsample score',
                        stat=stat,
                        phenod=pconfig.statusd(),
                        palette=pconfig.palette_diesase_status(),
                        figtype=figtype,
                        title=pconfig.phed([phe], kind='capital2')[phe],
                    )
                    fig.tight_layout(h_pad=-2, w_pad=-3)
                    # fig.tight_layout(h_pad=-2,w_pad=-5)
                    # fig.tight_layout(h_pad=-10,w_pad=-10)
                    # minus is allowed!!
                    # fig.tight_layout(h_pad=-0.05)
                    # fig.tight_layout(h_pad=-0.05)
                    # fig.tight_layout(h_pad=-1.0)
                    # fig.subplots_adjust(hspace=0.05)

                    # add dataset, vs
                    fname = droot / ('sample-score.summary.' + x_method + '-vs-' + y_method + '.' + dataset + '.' + stat + '.paper.figtype' + figtype)
                    logger.debug('fname {}'.format(fname))
                    pio.savefig_png_svg(fig, fname)
                    plt.close()
