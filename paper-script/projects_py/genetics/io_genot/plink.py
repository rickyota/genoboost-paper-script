

from logging import getLogger
import numpy as np
import pandas as pd


logger = getLogger(__name__)


def load_sample(fgenot, format, raw=False):
    # raw: return original encoding ex. sex (1,2) not (0,1)
    chrom_split = '%' in fgenot
    if chrom_split:
        fgenot = fgenot.replace('%', '1')

    if format == 'plink2' or format == 'plink2vzs':
        fsample = fgenot + '.psam'
        sample = load_sample_plink2(fsample, raw)
    elif format == 'plink':
        raise NotImplementedError
    else:
        raise RuntimeError
    return sample


def load_sample_plink2(fsample, raw=False):
    sample = pd.read_csv(fsample, delim_whitespace=True, dtype=str)

    sample.columns = sample.columns.str.replace('#', '')
    col_map = {'IID': 'id', 'SEX': 'sex'}
    sample = sample.rename(columns=col_map)

    if not raw:
        # change encoding from
        # {male: 1, female: 2}
        # to {male: 1, female: 0} (ukb way)
        assert set(sample['sex'].unique()) == {'1', '2'}, 'Unspecified sex in sample: {}'.format(fsample)
        sample.loc[:, 'sex'] = (2 - sample.loc[:, 'sex'].astype(np.int32)).astype(str)
    return sample

# bad example
# This is hard to test with StringIO
# def load_sample_plink2(fgenot, chrom_split):
#    if chrom_split:
#        fgenot = fgenot.replace('%', '1')
#    fsample = fgenot + '.psam'
#    sample = pd.read_csv(fsample, delim_whitespace=True)
#
#    sample.columns = sample.columns.str.replace('#', '')
#    col_map = {'IID': 'id', 'SEX': 'sex'}
#    sample = sample.rename(columns=col_map)
#    return sample


# def load_snv(fgenot, format, chrom=None, compress=False):
def load_snv(fgenot, format, chrom=None):
    # FIXME: when fgenot is pathlib.Path
    # FIXME: fgenot should be Path
    chrom_split = '%' in fgenot
    if chrom_split:
        if chrom is None:
            raise RuntimeError('indicate chrom')
        fgenot = fgenot.replace('%', str(chrom))

    if format == 'plink2vzs':
        fsnv = fgenot + '.pvar.zst'
        snv = load_snv_plink2(fsnv)
    elif format == 'plink2':
        fsnv = fgenot + '.pvar'
        snv = load_snv_plink2(fsnv)
    elif format == 'plink':
        fsnv = fgenot + '.bim'
        snv = load_snv_plink(fsnv)
    else:
        raise RuntimeError('Unknown format: ', format)
    return snv


def snv_col_map_plink2():
    col_map = {'CHROM': 'chrom', 'POS': 'pos', 'ID': 'id', 'REF': 'ref', 'ALT': 'alt',
               'ALT_FREQS': 'alt_frq', 'OBS_CT': 'obs_ct'  # .afreq
               }
    return col_map


# usually header=True
def load_snv_plink2(fsnv, header=True):
    # read_csv can read compressed .zst
    if header:
        snv = pd.read_csv(fsnv, header=0, delim_whitespace=True, dtype=str)
        snv.columns = snv.columns.str.replace('#', '')
        col_map = snv_col_map_plink2()
        # col_map = {'CHROM': 'chrom', 'POS': 'pos', 'ID': 'id', 'REF': 'ref', 'ALT': 'alt'}
        snv = snv.rename(columns=col_map)
    else:
        snv = pd.read_csv(fsnv, header=None,
                          names='chrom,pos,id,ref,alt'.split(','),
                          delim_whitespace=True, dtype=str)
    return snv


def load_snv_plink(fsnv):
    snv = pd.read_csv(fsnv, header=None,
                      names='chrom,id,gpos,pos,ref,alt'.split(','),
                      delim_whitespace=True, dtype=str)
    return snv


def load_freq_plink2(ffreq, header=True):
    if header:
        freq = pd.read_csv(ffreq, header=0, delim_whitespace=True, dtype=str)
        freq.columns = freq.columns.str.replace('#', '')
        col_map = snv_col_map_plink2()
        # col_map = {'CHROM': 'chrom', 'POS': 'pos', 'ID': 'id', 'REF': 'ref', 'ALT': 'alt', 'ALT_FREQS': 'alt_freq', 'OBS_CT': 'obs_ct'}
        freq = freq.rename(columns=col_map)
    else:
        # ex. after extract lines and no header
        freq = pd.read_csv(ffreq, header=None,
                           names='chrom,id,ref,alt,alt_frq,obs_ct'.split(','),
                           delim_whitespace=True, dtype=str)

    # print('freq',freq.info())
    freq = freq.astype({'alt_frq': np.float64, 'obs_ct': np.int64})
    # print('freq',freq.info())
    # this cannot overwite dtype
    # freq.loc[:, 'alt_frq'] = freq['alt_frq'].astype(np.float64)
    # freq.loc[:, 'obs_ct'] = freq['obs_ct'].astype(np.int64)

    return freq


# FIXME: snplist should have ('id', 'ref', 'alt')  or ('id','A1','A2')
def load_genot_rust(fgenot, genot_format, snplist):
    snv1 = load_snv(fgenot, genot_format)
