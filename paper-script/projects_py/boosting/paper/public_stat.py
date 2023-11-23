
from . import config as paper_config

def stat_to_paper_public(stat, cols):
    stat_public = stat.copy()
    stat_public = stat_public.loc[:, cols]

    if 'method' in cols:
        d= paper_config._methodd_capital()
        stat_public.loc[:, 'method'] = stat_public['method'].apply(lambda x:d[x])
    if 'phe' in cols:
        d=paper_config._phed_capital()
        # use original for simu ('random_0')
        stat_public.loc[:, 'phe'] = stat_public['phe'].apply(lambda x:d.get(x,x))
        #stat_public.loc[:, 'phe'] = stat_public['phe'].apply(lambda x:d[x])
    if 'stat_on' in cols:
        d= paper_config.accd()
        d|={'or_prop-1.0':'odds ratio at 1%',
            'or_prop-3.0':'odds ratio at 3%',
            'or_prop-5.0':'odds ratio at 5%',
            'or_prop-10.0':'odds ratio at 10%',
            }
        stat_public.loc[:, 'stat_on'] = stat_public['stat_on'].apply(lambda x:d[x])

    cold={'phe':'phenotype',
          'cvi': 'cross validation',
          'nsnv_use':'number of snvs',
          'stat_on':'metric'}
    stat_public=stat_public.rename(columns=cold)
        
    return stat_public
