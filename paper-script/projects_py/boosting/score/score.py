# move to filename_score.py ?

from logging import getLogger
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pathlib


# from ...genetics.io_genot import sample as sampleio


logger = getLogger(__name__)


# in sampleop
# def extract_score(score,dataset, cvi, cvs, sex='female'):
#    if sex == 'female':
#        pass
#    else:
#        pass


def parse_para(f):
    # parse para from file name
    # assume suffix only has one '.'
    # ok: '.score_withcov'
    # ng: '.score.withcov'
    # -> ex. ldpred_rho-0.1.wgt
    # -> now, allow .score.gz

    fpath = pathlib.Path(f)

    while fpath.suffix in {'.gz'}:
        fpath = fpath.with_suffix('')

    # base is str not Path
    base = fpath.stem
    pnames = base.split('_')

    parad = {}
    for pname in pnames:
        if '-' not in pname:
            continue
        ps = pname.split('-')

        parad[ps[0]] = ps[1]
    return parad


# def load_score(f):
#    dtype = defaultdict(lambda: np.float64)
#    dtype['id'] = str
#    scorecov = pd.read_csv(f, sep='\t', dtype=dtype)
#    return scorecov
# def load_score_phe(fscore, fphe, phe_name):
#    score = load_score(fscore)
#    phe = sampleio.load_phe(fphe, phe_name)
#    score_phe = pd.merge(score, phe, on='id', how='left')
#    return score_phe


def get_score_cols(score):
    # or .startswith('score')
    return [x for x in score.columns if (x != 'id') and (x != 'status')]
    # return [x for x in score.columns if (x != 'id') and (x != 'phe')]


def get_n_para(col):
    if col == 'score':
        # return None
        # use -1 only for pd => could forget
        return -1
    if col.startswith('score_n-'):
        return int(col.replace('score_n-', ''))
    raise RuntimeError('Cannot parse n_para: ', col)
