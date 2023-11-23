
from logging import getLogger
import numpy as np
import pandas as pd
from collections import defaultdict

from ..sample import sampleio as sampleio

logger = getLogger(__name__)


def load_score(f):
    # str
    # if not f.endswith('.gz'):
    # pahtlib.Path
    if f.suffix != '.gz':
        raise RuntimeError('sample name should end with gz {}'.format(f))
    dtype = defaultdict(np.float64, id='str')
    score = pd.read_csv(f, sep='\t', dtype=dtype, compression='gzip')
    #dtype ={'id':str}
    #score = pd.read_csv(f, sep='\t', compression='gzip')
    #score = score.astype({'id': str})
    #print('score ',score)
    return score


def load_score_fpheno(fscore, fpheno, phe):
    score = load_score(fscore)
    pheno = sampleio.load_pheno_phe(fpheno, phe)

    score_phe = pd.merge(score, pheno, on='id', how='left')
    return score_phe


def load_score_pheno(fscore, pheno, phe):
    score = load_score(fscore)
    pheno = sampleio.extract_phe(pheno, phe)

    score_pheno = pd.merge(score, pheno, on='id', how='left')
    return score_pheno


def write_score(f, score):
    #if not fout.endswith('.gz'):
    if f.suffix != '.gz':
        raise RuntimeError('sample name should end with gz {}'.format(f))
    score.to_csv(f, sep='\t', index=False, compression='gzip')
