
from logging import getLogger

logger = getLogger(__name__)


# TODO: rename cvs->cv
def extract_dataset(pheno, dataset, cvi, cvs, sex):
    if sex == 'both':
        id = cvs.loc[cvs['cv' + str(cvi)] == dataset, :].set_index('id').index
    elif sex == 'female':
        id = cvs.loc[(cvs['cv' + str(cvi)] == dataset) & (cvs['sex'] == '0'), :].set_index('id').index
    else:
        raise RuntimeError('Unknown sex', sex)

    # could happen
    #if len(id) == 0:
        #raise RuntimeError('No samples in cv')

    # stay order
    # phe.index.name = 'index'
    pheno = pheno.reset_index(names='index', drop=False)
    pheno = pheno.set_index('id')
    pheno_dataset = pheno.loc[id, :]
    pheno_dataset = pheno_dataset.reset_index(drop=False)
    pheno_dataset = pheno_dataset.set_index('index')
    pheno_dataset = pheno_dataset.sort_index()

    pheno_dataset = pheno_dataset.reset_index(drop=True)
    return pheno_dataset.copy()


def extract_align(covar, pheno):
    covar_index = covar.set_index('id')

    covar_align = covar_index.loc[pheno.set_index('id').index, :]
    covar_align = covar_align.reset_index(drop=False)

    if not (covar_align['id'] == pheno['id']).all():
        raise NotImplementedError('Some index in phe is not in cov')

    return covar_align.copy()
