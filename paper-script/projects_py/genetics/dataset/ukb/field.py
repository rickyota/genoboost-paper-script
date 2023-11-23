import pandas as pd
import itertools
from logging import getLogger

logger = getLogger(__name__)


def get_cov_fieldd_all():
    fields = {
        'sex': ['31'],
        'age': ['34', '52'],
        'PC': ['22009'],
    }
    return fields


def get_cov_fields(cov_names):
    fieldd_all = get_cov_fieldd_all()

    for cov_name in cov_names:
        if cov_name not in fieldd_all:
            RuntimeError('Unknown cov_name: ', cov_name)

    fieldd = {k: fieldd_all[k] for k in cov_names}
    return fieldd


def cov_fields_list(cov_names):
    fieldd = get_cov_fields(cov_names)
    fields_nest = list(fieldd.values())
    fields = list(itertools.chain.from_iterable(fields_nest))
    return fields


def load_ukb_cov(fukb, cov_names):
    fields_use = cov_fields_list(cov_names)
    ukb = load_ukb(fukb, fields_use)
    return ukb


def load_ukb(fukb, fields=[]):
    logger.debug('fields\n{}'.format(str(fields)))

    #def usecols(x):
    #    return (x == 'eid') or any([x.startswith(df) for df in fields])
    def usecols(x):
        return (x == 'eid') or any([x.startswith(df + '-') for df in fields])
    df = pd.read_csv(fukb,
                     sep=',',
                     header=0,
                     usecols=usecols,
                     dtype=str
                     )

    df = df.rename(columns={'eid': 'id'})

    return df
