import time
from logging import getLogger

logger = getLogger(__name__)


def add_sid(df, chrom='chrom', pos='pos'):
    df = df.assign(sid=df.apply(lambda x: '%s:%s' % (x[chrom], x[pos]), axis=1))
    # df = df.assign(sid=df.apply(lambda x: '%s:%s' % (x['chrom'], x['pos']), axis=1))
    return df


def add_sida(df, chrom='chrom', pos='pos', ref='a2', alt='a1'):
    # (chrom):(pos):(alt):(ref)
    df = df.assign(sida=df.apply(lambda x: '%s:%s:%s:%s' % (x[chrom], x[pos], x[alt], x[ref]), axis=1))
    # df = df.assign(sida=df.apply(lambda x: '%s:%s:%s:%s' % (x['chrom'], x['pos'], x['a1'], x['a2']), axis=1))
    return df


def add_sidannot(df, chrom='chrom', pos='pos', ref='a2', alt='a1'):
    # (chrom):(pos):(ref)>(alt)
    df = df.assign(sida=df.apply(lambda x: '%s:%s:%s>%s' % (x[chrom], x[pos], x[ref], x[alt]), axis=1))
    # df = df.assign(sida=df.apply(lambda x: '%s:%s:%s:%s' % (x['chrom'], x['pos'], x['a1'], x['a2']), axis=1))
    return df


def add_sids(df, chrom='chrom', pos='pos', ref='a2', alt='a1'):
    # https://stackoverflow.com/questions/11858472/string-concatenation-of-two-pandas-columns/54298586#54298586

    df = add_sid(df, chrom, pos)
    df = add_sida(df, chrom, pos, ref, alt)
    # df = add_sidannot(df,chrom,pos,alt,ref)

    # too much memory error
    # c0 = np.char.array(df['chrom'].to_numpy(dtype=str))
    # c1 = np.char.array(df['pos'].to_numpy(dtype=str))
    # c2 = np.char.array(df['A1'].to_numpy(dtype=str))
    # c3 = np.char.array(df['A2'].to_numpy(dtype=str))
    # df.loc[:, 'sida'] = (c0 + ':' + c1 + ':' + c2 + ':' + c3).astype(str)

    # print('sida', df.head())

    # very slow
    # + should be fast, so astype is slow?
    # df.loc[:, 'sida'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['A1'] + ':' + df['A2']

    return df


def match_allele(snv1, snv2, both=False,
                 return_filter=False,
                 allow_flip=True,  # ref <-> alt
                 allow_reverse=True,  # 'A' <-> 'T', 'C' <-> 'G'
                 order=False,  # order='left' means snv2's order is aligned to snv1
                 ignore_index=False,  # returned index is ranged
                 col1_ref='ref', col1_alt='alt',
                 col2_ref='ref', col2_alt='alt',
                 verbose=False):
    """
    return extracted snvs and bool array of snvs
    the order of snvs stays when order=False
    snv2's order is aligned to snv1 when order='left'
    if both: return (bool array of snv1, bool array of snv2)

    True if the snp in snv1 is also in snv2,
    allowing flipped allele, reversed strand.

    Do not assume input snvs are sorted.
    """

    if not ((order is False) or (order == 'left')):
        raise RuntimeError('Unknown order', order)

    assert set(snv1.columns) >= {'chrom', 'pos', col1_ref, col1_alt}
    assert set(snv2.columns) >= {'chrom', 'pos', col2_ref, col2_alt}
    # assert set(snv1.columns) >= {'chrom', 'pos', 'ref', 'alt'}
    # assert set(snv2.columns) >= {'chrom', 'pos', 'ref', 'alt'}

    # save original cols
    snv1_cols_ori = snv1.columns
    snv2_cols_ori = snv2.columns

    # assume temp col to use in this function is not in original snvs
    cols_tmp_use = ['use', 'match_pos', 'refr', 'altr', 'index_align']
    assert set(snv1.columns).isdisjoint(set(cols_tmp_use))
    assert set(snv2.columns).isdisjoint(set(cols_tmp_use))

    t0 = time.time()

    snv1 = snv1.copy()
    snv2 = snv2.copy()

    if verbose:
        logger.debug('snv1\n{}'.format(snv1.head()))
        logger.debug('snv2\n{}'.format(snv2.head()))

    # use ranged index
    snv1 = snv1.reset_index(drop=True, names='index').reset_index(drop=False)
    snv2 = snv2.reset_index(drop=True, names='index').reset_index(drop=False)

    # exclude snv not match to shorten time
    # judged by chrom and pos
    snv1.loc[:, 'match_pos'] = False
    snv2.loc[:, 'match_pos'] = False

    # set chrom and pos as index
    snv1 = snv1.set_index(['chrom', 'pos']).sort_index()
    snv2 = snv2.set_index(['chrom', 'pos']).sort_index()

    index_inter = snv1.index.unique().intersection(snv2.index.unique())
    snv1.loc[index_inter, 'match_pos'] = True
    snv2.loc[index_inter, 'match_pos'] = True

    if allow_reverse:
        nt_revd = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        snv1.loc[:, 'refr'] = snv1[col1_ref].map(nt_revd)
        snv1.loc[:, 'altr'] = snv1[col1_alt].map(nt_revd)
        # snv1.loc[:, 'refr'] = snv1['ref'].map(nt_revd)
        # snv1.loc[:, 'altr'] = snv1['alt'].map(nt_revd)

    snv1.loc[:, 'use'] = False
    snv2.loc[:, 'use'] = False

    # list all candidates of sida of wgt and match to bim
    snv2 = snv2.reset_index(drop=False).set_index(['chrom', 'pos', col2_ref, col2_alt]).sort_index()
    # snv2 = snv2.reset_index(drop=False).set_index(['chrom', 'pos', 'ref', 'alt']).sort_index()
    # for col_a1, col_a2 in (('A1', 'A2'), ('A1r', 'A2r'),  # match, reversed match
    # ('A2', 'A1'), ('A2r', 'A1r')):  # flipped match, flipped reversed match

    # col1 for snv1
    if allow_flip and allow_reverse:
        cols_pair = ((col1_ref, col1_alt), ('refr', 'altr'),  # match, reversed match
                     (col1_alt, col1_ref), ('altr', 'refr'))  # flipped match, flipped reversed match
        # cols_pair = (('ref', 'alt'), ('refr', 'altr'),  # match, reversed match
        # ('alt', 'ref'), ('altr', 'refr'))  # flipped match, flipped reversed match
    elif allow_flip and (not allow_reverse):
        cols_pair = ((col1_ref, col1_alt),   # match
                     (col1_alt, col1_ref))  # flipped match
        # cols_pair = (('ref', 'alt'),   # match
        # ('alt', 'ref'))  # flipped match
    elif (not allow_flip) and allow_reverse:
        raise RuntimeError('Usually not applicable')
        cols_pair = (('ref', 'alt'), ('refr', 'altr'))  # match, reversed match
    elif (not allow_flip) and (not allow_reverse):
        # cols_pair = (('ref', 'alt'))  # match
        cols_pair = ((col1_ref, col1_alt))  # match

    if order == 'left':
        snv2.loc[:, 'index_align'] = -1

    for col_a1, col_a2 in cols_pair:
        if verbose:
            logger.debug("col_a1, col_a2={}, {}".format(col_a1, col_a2))
            # print("wgt", snv1.info())

        snv1 = snv1.reset_index(drop=False).set_index(['chrom', 'pos', col_a1, col_a2]).sort_index()

        # wgt['sida_'] = wgt['chrom'].apply(
        # str) + ':' + wgt['pos'].apply(str) #+ ':' + wgt[col_a1].apply(str) + ':' + wgt[col_a2].apply(str)
        # wgt.loc[:,'sida_'] =  ':' + wgt[col_a1] + ':' + wgt[col_a2]
        # snv1.loc[:, 'sida_'] = snv1['chrom'].astype(
        #    str) + ':' + snv1['pos'].astype(str) + ':' + snv1[col_a1] + ':' + snv1[col_a2]

        # snv1 = snv1.set_index('sida_')

        if verbose:
            logger.debug("wgt {}".format(snv1.head()))
            logger.debug("bim {}".format(snv2.head()))

        # index_inter_ = snv1.index.intersection(snv2.index)
        index_inter = snv1[snv1['match_pos']].index.intersection(snv2[snv2['match_pos']].index)

        logger.debug("intersect number {}".format(len(index_inter)))

        if verbose:
            pass
            # print("intersection", index_inter)

        snv1.loc[index_inter, 'use'] = True
        snv2.loc[index_inter, 'use'] = True

        if order == 'left':
            snv2.loc[index_inter, 'index_align'] = snv1.loc[index_inter, 'index']

        if verbose:
            pass
            # print("sida_ok", sida_ok)
            # print("sida_ok len", len(sida_ok))

        # snv1 = snv1.reset_index(drop=True)

    snv1 = snv1.reset_index(drop=False).set_index('index').sort_index().reset_index(drop=True)
    snv2 = snv2.reset_index(drop=False).set_index('index').sort_index().reset_index(drop=True)

    filter1 = snv1['use'].to_numpy()
    filter2 = snv2['use'].to_numpy()

    snv1 = snv1.loc[filter1, :]
    # snv1 = snv1.loc[filter1, snv1_cols_ori]
    snv2 = snv2.loc[filter2, :]
    # snv2 = snv2.loc[filter2, snv2_cols_ori]
    if order == 'left':
        if (snv2['index_align'] == -1).any():
            raise RuntimeError('Sth wrong')
        snv2 = snv2.sort_values('index_align')

    if ignore_index:
        snv1 = snv1.reset_index(drop=True)
        snv2 = snv2.reset_index(drop=True)

    snv1 = snv1.loc[:, snv1_cols_ori]
    snv2 = snv2.loc[:, snv2_cols_ori]

    logger.debug("snv1 intersect {}".format(len(snv1)))
    logger.debug("snv2 intersect {}".format(len(snv2)))

    logger.debug('time of isin_genot_sida {}'.format(time.time() - t0))

    assert len(snv1) == len(snv2)

    if both:
        if return_filter:
            return snv1, filter1, snv2, filter2
        else:
            return snv1, snv2
    else:
        if return_filter:
            return snv1, filter1
        else:
            return snv1
