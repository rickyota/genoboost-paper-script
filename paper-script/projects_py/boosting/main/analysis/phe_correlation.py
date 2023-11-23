
import os
import argparse
import pathlib
import pandas as pd
import numpy as np
from sklearn import metrics
from scipy import stats as scistats
# from scipy.stats import contingency
from logging import getLogger

# from ....genetics.io_genot import plink
# from ....genetics.snv import snv as snvop
# from ....genetics.snv import snvio
from ....genetics.sample import sampleio as sampleio
from ....genetics.sample import sample as sampleop
# from ...genetics.io_genot import sample as sampleio
# from ...genetics.io_genot import score as scoreio
# from ...genetics.sample import sample as sampleop
# from ..io import filename as iof
# from ..score import score as scoreop
# from ...wgt import wgt as wgtop
# from ..validate.paper import auc as validate_paper_auc
from .plot import phes_correlation as plot_phe_correlation

from ....system.logger import logger_setting

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--dout',
                        help='output dir')

    # ./result/nonaddpgs/dataset.12/
    parser.add_argument('--dresult',
                        help='output dir')

    parser.add_argument('--ddata',
                        help='output dir')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--phes',
                        help='')

    parser.add_argument('--phes-female',
                        help='')

    parser.add_argument('--fcv',
                        help='')

    parser.add_argument('--cvi',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def get_phes_bothsex(phes, phes_female):
    return [x for x in phes if x not in phes_female]


def get_phes_female(phes, phes_female):
    return [x for x in phes if x in phes_female]


def jaccard_similarity(phes1, phes2):
    stat = metrics.jaccard_score(phes1, phes2)
    return stat


def chisq_test(phes1, phes2):
    res = scistats.contingency.crosstab(phes1, phes2)
    count = pd.DataFrame(res.count, index=res.elements[0], columns=res.elements[1])
    # format order
    label_order = [0, 1]
    count = count.loc[label_order, label_order]
    # logger.debug('count: {}'.format(count))

    stat = scistats.chi2_contingency(count)
    stat = stat.pvalue
    # logger.debug('stat: {:e}'.format(stat))

    # do this in plot
    #if stat==0.0:
        # cannot use np.float128, so should use smallest float
        # TMP
        #stat= np.nextafter(0,1)

    return stat


def phes_cor(pheno, stat_fn, phes, phes_female, cvs, cvi, sex):
    if sex == 'both':
        phes_use = get_phes_bothsex(phes, phes_female)
    elif sex == 'female':
        phes_use = get_phes_female(phes, phes_female)
    else:
        raise NotImplementedError

    pheno = sampleop.extract_dataset(pheno, dataset='ts', cvi=cvi, cvs=cvs, sex=sex)

    pheno_use = pheno.loc[:, phes_use].copy()
    # pheno_use = pheno.loc[samples, phes_use].copy()

    stats = pd.DataFrame([], index=phes_use, columns=phes_use, dtype=float)
    for phe1 in phes_use:
        for phe2 in phes_use:
            stat = stat_fn(pheno_use[phe1], pheno_use[phe2])
            # stat = jaccard_similarity(pheno_use[phe1], pheno_use[phe2])

            # TODO: chi-square test p-value
            stats.loc[phe1, phe2] = stat

    return stats


def phe_cor(dout, dresult, ddata, fpheno, phes, phes_female, fcv, cvi):
    # correlation of phes_both for test samples for both sex

    cvs = sampleio.load_cv(fcv)
    pheno = sampleio.load_pheno(fpheno)

    sex = 'both'
    droot = dout / 'phe_correlation' / sex
    os.makedirs(droot, exist_ok=True)

    stat_fns = [('jaccard_similarity', jaccard_similarity),
                ('chisq', chisq_test)]

    for (stat_name, stat_fn) in stat_fns:
        stats = phes_cor(pheno, stat_fn, phes, phes_female, cvs, cvi, sex)
        # logger.debug('stat: {},{}'.format(stat_name, stats.head()))
        logger.debug('stat: {},{}'.format(stat_name, stats.to_string()))
        plot_phe_correlation.plot_phe_correlation(droot, stats, stat_name)
        # plot_phe_correlation.plot_phe_correlation(droot, stats, 'jaccard_similarity')


def main():
    args = argument()

    phes = args.phes.split(' ')
    logger.debug('phes: {}'.format(phes))
    phes_female = args.phes_female.split(' ')

    phe_cor(pathlib.Path(args.dout),
            pathlib.Path(args.dresult),
            pathlib.Path(args.ddata),
            pathlib.Path(args.fpheno),
            phes,
            phes_female,
            pathlib.Path(args.fcv),
            int(args.cvi),
            )


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
