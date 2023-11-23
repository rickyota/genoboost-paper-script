
from logging import getLogger
import pandas as pd

logger = getLogger(__name__)


def write_sample(fout, sample):
    """
    use 'IID' not 'id' so that the fsample fits the plink input format
     """
    sample = sample.copy()

    col_map = {'id': 'IID'}
    # then use plink.write_sample()
    # col_map = {'id': 'IID','sex': 'SEX'}
    sample = sample.rename(columns=col_map)
    sample.to_csv(fout, sep='\t', index=False)


def load_cv(fin):
    df = pd.read_csv(fin, header=0, delim_whitespace=True, dtype=str)
    col_map_in = {'IID': 'id'}
    df = df.rename(columns=col_map_in)
    return df


def load_pheno(fin):
    df = pd.read_csv(fin, header=0, delim_whitespace=True)
    col_map_in = {'IID': 'id'}
    df = df.rename(columns=col_map_in)
    df = df.astype({'id': str})
    return df


def extract_phe(pheno, phe):
    # print('pheno', pheno)
    pheno = pheno[['id', phe]].copy()
    # rename phe_name -> 'status'
    # use 'status' to avoid confusion with disease name
    col_map_in = {phe: 'status'}
    # rename phe_name -> 'phe'
    #col_map_in = {phe: 'phe'}
    pheno = pheno.rename(columns=col_map_in)
    return pheno


def load_pheno_phe(fin, phe):
    pheno_whole = load_pheno(fin)
    pheno = extract_phe(pheno_whole, phe)
    return pheno


# usually, id+covs+phes
# def load_pheno_whole_str(fin):
#    df = pd.read_csv(fin, header=0, delim_whitespace=True, dtype=str)
#    col_map_in = {'IID': 'id'}
#    df = df.rename(columns=col_map_in)
#    return df
#
#
# def extract_phe(pheno, phe):
#    print('pheno', pheno)
#    pheno = pheno[['id', phe]].copy()
#    # rename phe_name -> 'phe' ?
#    col_map_in = {phe: 'phe'}
#    pheno = pheno.rename(columns=col_map_in)
#    pheno = pheno.astype({'phe': np.int64})
#    return pheno
#
#
# def load_pheno_phe(fin, phe):
#    pheno_whole = load_pheno_whole_str(fin)
#    pheno = pheno_whole[['id', phe]].copy()
#    # rename phe_name -> 'phe' ?
#    col_map_in = {phe: 'phe'}
#    pheno = pheno.rename(columns=col_map_in)
#    pheno = pheno.astype({'phe': np.int64})
#    return pheno
#

def parse_cov_names(cov_columns, cols):
    # cov_columns: 'age,PC_1-PC_10'
    # return ['age','PC_1,'PC_2',...,'PC_10']

    if len(cols) != len(set(cols)):
        raise RuntimeError('Same col is in cols.')

    cov_names = []
    for cov_col in cov_columns.split(','):
        if '-' in cov_col:
            xs = cov_col.split('-')
            assert len(xs) == 2
            x0, x1 = xs[0], xs[1]
            assert '_' in x0
            assert '_' in x1

            index0 = cols.index(x0)
            index1 = cols.index(x1)
            cov_name = cols[index0:(index1 + 1)]
            cov_names += cov_name
        else:
            cov_name = cov_col
            cov_names += [cov_name]
    return cov_names


def load_cov(fpheno,
             cov_names=None,
             cov_names_unparse=None):
    if cov_names is None and cov_names_unparse is None:
        raise RuntimeError('Use either cov_names or cov_names_unparse')

    if cov_names is not None and not isinstance(cov_names, list):
        raise RuntimeError('cov_names should be list')

    pheno_whole = load_pheno(fpheno)

    if cov_names_unparse is not None:
        cov_names = parse_cov_names(cov_names_unparse, list(pheno_whole.columns))

    logger.debug('cov_names {}'.format(cov_names))

    covar = pheno_whole[['id'] + cov_names].copy()
    # create new one to change dtype
    # covar_new = covar.loc[:, ['id']].copy()
    # for col in cov_names:
    #    # print(pd.to_numeric(cov[col]).info())
    #    covar_new.loc[:, col] = pd.to_numeric(covar[col])
    # logger.debug('cov_new\n{}'.format(covar_new.info()))
    logger.debug('cov info: {}'.format(covar.info()))
    #logger.debug('cov: {}'.format(covar.info()))
    return covar

# use str
# def load_cov(fpheno,
#             cov_names=None,
#             cov_names_unparse=None):
#    if cov_names is None and cov_names_unparse is None:
#        raise RuntimeError('Use either cov_names or cov_names_unparse')
#
#    if cov_names is not None and not isinstance(cov_names, list):
#        raise RuntimeError('cov_names should be list')
#
#    pheno_whole = load_pheno_whole_str(fpheno)
#
#    if cov_names_unparse is not None:
#        cov_names = parse_cov_names(cov_names_unparse, list(pheno_whole.columns))
#
#    logger.debug('cov_names {}'.format(cov_names))
#
#    covar = pheno_whole[['id'] + cov_names].copy()
#    # create new one to change dtype
#    covar_new = covar.loc[:, ['id']].copy()
#    for col in cov_names:
#        # print(pd.to_numeric(cov[col]).info())
#        covar_new.loc[:, col] = pd.to_numeric(covar[col])
#    # ng
#    # cov_new.loc[:, cov_names] = cov[cov_names].convert_dtypes()
#    # cov_new.loc[:, cov_names] = cov[cov_names].convert_dtypes(convert_string=False,convert_boolean=False)
#    logger.debug('cov_new\n{}'.format(covar_new.info()))
#    return covar_new
