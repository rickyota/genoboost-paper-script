
import argparse
import os
import pathlib
import pandas as pd
import numpy as np
from logging import getLogger

from ....genetics.io_genot import plink
from ....genetics.snv import snv as snvop
from ....genetics.snv import snvio
# from ...genetics.io_genot import sample as sampleio
# from ...genetics.io_genot import score as scoreio
# from ...genetics.sample import sample as sampleop
# from ..io import filename as iof
# from ..score import score as scoreop
from ...wgt import wgt as wgtop
from ...ss import ss as ssop
# from ..validate.paper import auc as validate_paper_auc

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

    # use --para-regon
    # parser.add_argument('--regon-ss',
    #                    help='')

    # parser.add_argument('--method',
    #                    help='')

    parser.add_argument('--mid-path',
                        nargs="+",
                        help='')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--phe',
                        help='')

    parser.add_argument('--cov-names',
                        # nargs='+',
                        help='ex. age,sex,PC_1-PC_10. When use "-", the cov ("PC") should be the same.')

    parser.add_argument('--fcv',
                        help='')

    parser.add_argument('--cvi',
                        help='')

    parser.add_argument('--fwgt',
                        help='')

    parser.add_argument('--fmap-hg38',
                        help='')

    parser.add_argument('--nsnv',
                        help='')

    parser.add_argument('--fgenot',
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2',
                        help='plink2vzs, plink2, plink, etc.')

    parser.add_argument('--ffreq',
                        help='')

    parser.add_argument('--fcatalog',
                        help='')

    parser.add_argument('--fgene',
                        help='')

    parser.add_argument('--fgene_rna',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


# TODO: mv somewhere: wgtop
def merge_freq(wgt, ffreq, use_sid_tmp=False):
    freq = plink.load_freq_plink2(ffreq)
    logger.debug('freq: {}'.format(freq.info()))
    logger.debug('freq: {}'.format(freq.head()))
    # freq_check = freq[['id', 'ref', 'alt']].set_index('id')
    # snv_check = wgt.set_index(('id'))
    # if (freq_check['ref'] != snv_check['a2']).any() or (freq_check['alt'] != snv_check['a1']).any():
    #    raise RuntimeError('unmatch allele')

    if use_sid_tmp:
        logger.debug('USE SID')
        freq = snvop.add_sid(freq)
        freq_check = freq[['sid', 'ref', 'alt']].copy()
        wgt_check = wgt[['sid', 'a1', 'a2']].copy()
        wgt_check = pd.merge(wgt_check, freq_check, on='sid', how='left')
        logger.debug('wgt_check: {}'.format(wgt_check.head()))

        if (wgt_check['a1'] != wgt_check['alt']).any() or (wgt_check['a2'] != wgt_check['ref']).any():
            logger.error('unmatch: {}'.format(wgt_check[wgt_check['a1'] != wgt_check['alt']]))
            raise RuntimeError('unmatch allele')

        freq = freq[['sid', 'alt_frq']].rename({'alt_frq': 'a1_frq'})
        wgt = pd.merge(wgt, freq, on='sid', how='left')

    else:
        freq_check = freq[['id', 'ref', 'alt']].copy()
        wgt_check = wgt[['id', 'a1', 'a2']].copy()
        wgt_check = pd.merge(wgt_check, freq_check, on='id', how='left')
        logger.error('wgt_check: {}'.format(wgt_check.head()))

        if (wgt_check['a1'] != wgt_check['alt']).any() or (wgt_check['a2'] != wgt_check['ref']).any():
            logger.error('unmatch: {}'.format(wgt_check[wgt_check['a1'] != wgt_check['alt']]))
            raise RuntimeError('unmatch allele')

        freq = freq[['id', 'alt_frq']].rename({'alt_frq': 'a1_frq'})
        wgt = pd.merge(wgt, freq, on='id', how='left')

    return wgt


def phed_catalog():
    # TODO: cad include 'coronary heart disease'??
    # lower cap
    phenotyped = {'t2d': 'type 2 diabetes', 'cad': 'coronary artery disease', 'bc': 'breast cancer', 'af': 'atrial fibrillation',
                  'ibd': 'inflammatory bowel disease',
                  'ra': 'rheumatoid arthritis',
                  'gout': 'gout',
                  'ad': "alzheimer's disease (late onset)",
                  'acd': "alzheimer's disease (late onset)",  # use ad instead
                  # 'acd':'dementia', ######## NO assoc study for all-cause dementia
                  'cc': 'colorectal cancer',
                  'atm': 'asthma',
                  'psr': 'psoriasis',
                  }
    return phenotyped


def load_catalog(fcatalog, phe):
    phed = phed_catalog()
    phe_catalog = phed[phe]

    catalog = pd.read_csv(fcatalog, sep='\t', header=0, low_memory=False, encoding='utf-8')

    # print("catalog", catalog)
    # print("catalog", catalog.info())

    catalog = catalog[catalog['DISEASE/TRAIT'].str.lower() == phe_catalog]
    # catalog = catalog[catalog['DISEASE/TRAIT'].str.lower().str.contains(phenotype_catalog)]

    if len(catalog) == 0:
        logger.error('phenotype not found: {}'.format(phe))
        raise RuntimeError('No catalog of disease: ', phe_catalog)

    catalog = catalog.rename(columns={'CHR_ID': 'chrom', 'CHR_POS': 'pos'})

    # exclude exception
    catalog = catalog.dropna(subset=['pos'])
    catalog = catalog[~catalog['pos'].str.contains('x')].copy()

    # pos might be '123;124;125;'
    def f_pos(x):
        if ';' not in x:
            return x
        else:
            return x.split(';')[0]
    catalog.loc[:, 'pos'] = catalog['pos'].map(f_pos)
    catalog = catalog.astype({'pos': int})
    # chrom=str
    # catalog = catalog.astype({'chrom': int})

    logger.debug("catalog_pos: {}".format(catalog.info()))

    logger.debug("catalog: {}".format(catalog[['FIRST AUTHOR', 'DISEASE/TRAIT', 'chrom', 'pos']]))

    if len(catalog) == 0:
        raise RuntimeError('No snv in catalog after qc.')

    return catalog


def map_wgt_to_catalog(wgt, catalog):

    wgt.loc[:, 'catalog'] = False
    wgt.loc[:, 'catalog_rel'] = 'NaN'
    wgt.loc[:, 'catalog_gene'] = 'NaN'

    # 1Mbp
    bp_catalog_nearest = 1000000

    for ind in wgt.index:
        # use hg38
        chrom, pos = wgt.loc[ind, ['chrom_hg38', 'pos_hg38']]
        # logger.debug("chrom,pos: {}, {}".format(chrom, pos))

        catalog_chrom = catalog[catalog['chrom'] == str(chrom)]
        # logger.debug("catalog: {}".format(catalog_chrom[['chrom','pos']]))

        if len(catalog_chrom) == 0:
            continue

        ind_nearest = catalog_chrom.index[(catalog_chrom['pos'] - pos).abs().argsort()][0]
        # logger.debug("ind_nearest: {}".format(ind_nearest))

        diff = pos - int(catalog_chrom.loc[ind_nearest, 'pos'])
        # print("pos",int(catalog_chrom.loc[ind_nearest,'pos']))
        # print("diff",diff)

        if np.abs(diff) < bp_catalog_nearest:
            wgt.loc[ind, 'catalog'] = True
            wgt.loc[ind, 'catalog_rel'] = str(diff) + 'bp' if diff <= 0 else '+' + str(diff) + 'bp'
            wgt.loc[ind, 'catalog_gene'] = str(catalog_chrom.loc[ind_nearest, 'MAPPED_GENE'])
            # wgt.loc[ind,'catalog_gene']=catalog_chrom.loc[ind_nearest,'MAPPED_GENE']

    return wgt


def annotate_gene_overlap(wgt, gene, is_rna=False):
    # should allow overlapping genes
    # if snv is in several genes, write all of them
    # if not in any gene, write nearest gene

    bp_gene_nearest = 100000

    gene_name = 'gene'
    gene_rel = 'gene_rel'
    gene_biotype = 'gene_biotype'
    if is_rna:
        gene_name += '_rna'
        gene_rel += '_rna'
        gene_biotype += '_rna'

    wgt.loc[:, gene_name] = 'NaN'
    wgt.loc[:, gene_rel] = 'NaN'
    if 'biotype' in gene:
        wgt.loc[:, gene_biotype] = 'NaN'

    # wgt.loc[:, 'gene'] = 'NaN'
    # wgt.loc[:, 'gene_rel'] = 'NaN'
    # if 'biotype' in gene:
    #    wgt.loc[:, 'gene_biotype'] = 'NaN'

    # regard +- 100kb as 'gene'
    # put 'in' in 'gene_rel' if in the gene
    # put '+123bp'

    for ind in wgt.index:
        chrom, pos = wgt.loc[ind, ['chrom', 'pos']]

        gene_chrom = gene[gene['chrom'] == str(chrom)]
        # gene_chrom = gene[gene['chrom'] == chrom]

        # if len(gene_chrom)==0

        # check if in gene
        gene_in = gene_chrom[(gene_chrom['start'] <= pos) & (pos < gene_chrom['end'])]
        if len(gene_in) > 0:
            # snv could be in several genes
            wgt.loc[ind, gene_name] = gene_in['gene'].str.cat(sep=',')
            wgt.loc[ind, gene_rel] = 'in'
            if 'biotype' in gene:
                wgt.loc[ind, gene_biotype] = gene_in['biotype'].str.cat(sep=',')
                # wgt.loc[ind, gene_biotype] = gene_in.iloc[0, :]['biotype']
            # wgt.loc[ind, 'gene'] = gene_in['gene'].str.cat(sep=',')
            # wgt.loc[ind, 'gene_rel'] = 'in'
            # if 'biotype' in gene:
                # wgt.loc[ind, 'gene_biotype'] = gene_in.iloc[0, :]['biotype']
        else:
            # search nearest
            # first start
            ind_nearest_start = gene_chrom.index[(gene_chrom['start'] - pos).abs().argsort()][0]
            ind_nearest_end = gene_chrom.index[(gene_chrom['end'] - pos).abs().argsort()][0]
            # ind_nearest_start=gene_chrom.index[(gene_chrom['start']-pos).abs().argsort()]
            # ind_nearest_start=ind_nearest_start[0]
            # ind_nearest_end=gene_chrom.index[(gene_chrom['end']-pos).abs().argsort()]
            # ind_nearest_end=ind_nearest_end[0]

            diff_start = pos - gene_chrom.loc[ind_nearest_start, 'start']
            diff_end = pos - gene_chrom.loc[ind_nearest_end, 'end']

            if np.abs(diff_start) >= bp_gene_nearest and np.abs(diff_end) >= bp_gene_nearest:
                continue

            if np.abs(diff_start) <= np.abs(diff_end):
                wgt.loc[ind, gene_name] = gene_chrom.loc[ind_nearest_start, 'gene']
                # wgt.loc[ind, 'gene'] = gene_chrom.loc[ind_nearest_start, 'gene']
                wgt.loc[ind, gene_rel] = str(diff_start) + 'bp'
                # wgt.loc[ind, 'gene_rel'] = str(diff_start) + 'bp'
                if 'biotype' in gene:
                    wgt.loc[ind, gene_biotype] = gene_chrom.loc[ind_nearest_start, 'biotype']
            else:
                wgt.loc[ind, gene_name] = gene_chrom.loc[ind_nearest_end, 'gene']
                wgt.loc[ind, gene_rel] = '+' + str(diff_end) + 'bp'
                if 'biotype' in gene:
                    wgt.loc[ind, gene_biotype] = gene_chrom.loc[ind_nearest_end, 'biotype']

    # print("wgt", wgt[['sida', 'gene', 'gene_rel']])
    logger.debug("wgt gene {}".format(wgt[['sida', gene_name, gene_rel]]))

    return wgt


def add_catalog(wgt, fcatalog, phe, fmap_hg38, use_sid_tmp=False):
    map_hg38 = snvio.load_fmap_hg38(fmap_hg38)
    # print(map_hg38)
    map_hg38 = snvop.add_sid(map_hg38, 'chrom_hg38', 'pos_hg38')
    map_hg38 = map_hg38.rename(columns={'sid': 'sid_hg38'})

    if use_sid_tmp:
        # sida in hg19
        map_hg38 = map_hg38.rename(columns={'id': 'sida'})
        wgt = wgt.merge(map_hg38[['chrom_hg38', 'pos_hg38', 'sid_hg38', 'sida']], on='sida', how='left')
        # wgt = wgt.merge(map_hg38[['chrom_hg38', 'pos_hg38', 'sid_hg38', 'id']], left_on='sida', right_on='id', how='left')
    else:
        # wgt = wgt.merge(map_hg38[['chrom_hg38', 'pos_hg38', 'sid_hg38']], on='id', how='left')
        wgt = wgt.merge(map_hg38[['chrom_hg38', 'pos_hg38', 'sid_hg38', 'id']], on='id', how='left')
#    logger.debug('wgt {}'.format(wgt.info()))

    catalog = load_catalog(fcatalog, phe)
    logger.debug('catalog {}'.format(catalog))
    logger.debug('catalog {}'.format(catalog.info()))

    logger.debug('map_wgt_to_catalog')
    wgt = map_wgt_to_catalog(wgt, catalog)

    # print(wgt)
    logger.debug('wgt {}'.format(wgt[['catalog', 'catalog_rel']]))

    return wgt


def add_rs(wgt, fgenot, genot_format):
    snv = plink.load_snv(fgenot, genot_format)
    snv = snvop.add_sid(snv)

    wgt = wgt.drop(columns=['id'])
    logger.debug('snv {}'.format(snv.info()))
    wgt = pd.merge(wgt, snv[['sid', 'id']], on='sid', how='left')
    return wgt


def calculate_odd(wgt, fgenot, genot_format, fpheno, phe, fcv, cov_names):
    # TODO
    pass


def _write_paper(fout, wgt):

    cold_show = {
        'boost_iteration_uniq': 'index',
        'boost_iteration': 'iteration',
        # 'boost_t_snv': 'iteration',
        'boost_kind': 'kind',
        'id': 'var',
        # 'vid': 'var',
        # 'boost_var_annot': 'variant',
        'chrom': 'chrom',
        'pos': 'pos',
        # TODO: sid_hg38
        'sid_hg38': 'variant_hg38',
        # 'boost_var_annot_hg38': 'variant_hg38',
        # 'boost_var':'variant', # 'variant_hg38
        'boost_score0': 's0_boost',
        'boost_score1': 's1_boost',
        'boost_score2': 's2_boost',
        'boost_eps': 'eps_boost',
        # 'score1': 'brscore1',
        # 'score2': 'brscore2',
        # 'boost_eps': 'effect_inf',
        # add or of sample count as well
        'alt_frq': 'MAF',
        'boost_model_est': 'model_boost',
        # 'model_est': 'bmodel',
        'model_est_aic': 'aic_model',
        # 'pval_ukb': 'additive p-value',
        # 'pval_ukb_rec': 'recessive p-value',
        # 'pval_ukb_dom': 'dominant p-value',
        # 'pval_ukb_het': 'hetonly p-value',

        'p_add': 'additive p-value',
        'p_rec': 'recessive p-value',
        'p_dom': 'dominant p-value',
        'p_het': 'hetonly p-value',

        'gene': 'gene',
        'gene_rel': 'gene_distance',
        'gene_rna': 'gene_rna',
        'gene_rel_rna': 'gene_distance_rna',
        'gene_biotype_rna': 'gene_biotype_rna',
    }

    cols_use = [x for x in cold_show if x in wgt.columns]
    wgt_paper = wgt[cols_use].rename(columns=cold_show)

    # format float for pval
    if 'additive p-value' in wgt:
        format_pval = '{:.2e}'
        wgt_paper.loc[:, 'additive p-value'] = wgt_paper['additive p-value'].map(lambda x: str(format_pval.format(x)))
        wgt_paper.loc[:, 'recessive p-value'] = wgt_paper['recessive p-value'].map(lambda x: str(format_pval.format(x)))
        wgt_paper.loc[:, 'dominant p-value'] = wgt_paper['dominant p-value'].map(lambda x: str(format_pval.format(x)))
        wgt_paper.loc[:, 'hetonly p-value'] = wgt_paper['hetonly p-value'].map(lambda x: str(format_pval.format(x)))

    wgt_paper.to_csv(fout, sep='\t', na_rep='NaN', index=False,
                     float_format='%.4f')


def _write_catalog(fout, wgt):

    cols_show = [
        'boost_iteration',
        'boost_iteration_uniq',
        'id', 'sid', 'sid_hg38',
        # 'boost_var_annot', 'boost_var_annot_hg38', 'vid',  # 'boost_var', #'pos_hg38',
        'catalog_rel', 'catalog_gene',
        'gene', 'gene_rel',
        'gene_rna', 'gene_rel_rna', 'gene_biotype_rna',
        'alt_frq',
        'score1', 'score2',
        'boost_model_est', 'model_est_aic',
        # 'model_est', 'model_est_aic',
        # 'pval_ukb', 'pval_ukb_rec', 'pval_ukb_dom', 'pval_ukb_het',
        'p_add', 'p_rec', 'p_dom', 'p_overdom',
        'pval_gwas',
        'boost_eps'
    ]
    cols_use = [x for x in cols_show if x in wgt.columns]
    # cols_use = [x for x in wgt.columns if x in cols_show]
    wgt_out = wgt[cols_use].copy()
    # TMP comment
    # cols = ['pval_ukb', 'pval_ukb_rec', 'pval_ukb_dom', 'pval_ukb_het', 'pval_gwas']
    # wgt_out.loc[:, cols] = wgt_out[cols].applymap(lambda x: '%.1e' % x)

    def f_gene(x):
        if len(x) < 30:
            return x
        else:
            return x[:30] + '...'

    # print(wgt_out.head())
    # print(wgt_out.info())
    wgt_out.loc[:, 'catalog_gene'] = wgt_out['catalog_gene'].map(f_gene)
    wgt_out.to_csv(fout,
                   sep='\t', na_rep='NaN', float_format='%.4f', index=False)


def write_paper(dout, wgt, prefix):

    os.makedirs(dout, exist_ok=True)

    _write_catalog(dout / ('catalog' + prefix + '.tsv'), wgt)
    _write_paper(dout / ('paper' + prefix + '.tsv'), wgt)

    # _write_catalog(dout / 'catalog.allsnv.tsv', wgt)
    # _write_paper(dout / 'paper.allsnv.tsv', wgt)


def add_gene(wgt, fgene, fgene_rna):
    # add gene
    gene = snvio.load_gene(fgene)
    logger.debug('annotate_gene_overlap')
    wgt = annotate_gene_overlap(wgt, gene)
    # logger.debug('wgt {}'.format(wgt))

    # add gene or rna
    # gene = snvio.load_gene(fgene)
    # logger.debug('gene {}'.format(gene))
    gene_rna = snvio.load_gene(fgene_rna, is_rna=True)
    logger.debug('gene_rna {}'.format(gene_rna))
    logger.debug('annotate_gene_overlap')
    wgt = annotate_gene_overlap(wgt, gene_rna, is_rna=True)
    # rename to _rna

    return wgt


def format_wgt(wgt):
    wgt.loc[:, 'boost_iteration_uniq'] = range(len(wgt))
    return wgt


def extract_nonadd_table(wgt, is_rna=False):
    def fn_nonadd_table(x):
        # True; if catalog_gene is None and model_est is not add and gene in not NaN and gene_rel is in or <10bp

        if is_rna:
            gene_col = 'gene_rna'
            gene_rel_col = 'gene_rel_rna'
        else:
            gene_col = 'gene'
            gene_rel_col = 'gene_rel'

        # not in catalog
        # if not np.isnan(x['catalog_gene']):
        # if not (x['catalog_gene'] == 'NaN'):
        # TODO: should be np.nan not 'NaN'
        if x['catalog_gene'] != 'NaN':
            return False

        # not add
        # if x['model_est'] == 'add':
        if x['boost_model_est'] == 'add':
            return False

        # if not np.isnan(x['gene']):
        # if x['gene'] == 'NaN':
        if x[gene_col] == 'NaN':
            return False

        def is_gene_rel_near(y):
            # return True if 'in' or <10kb
            # if np.isnan(y):
            if y == 'NaN':
                return False
            if y == 'in':
                return True
            if not y.endswith('bp'):
                raise RuntimeError('why not ends with bp?')

            y_bp = int(y[:-2])

            # WITHIN 10 Kbp
            if np.abs(y_bp) < 10000:
                return True
            else:
                return False

        # if not is_gene_rel_near(x['gene_rel']):
        if not is_gene_rel_near(x[gene_rel_col]):
            return False

        return True

    filt = wgt.apply(fn_nonadd_table, axis=1)
    return wgt.loc[filt, :]


def extract_wgt_first_in_ld_width(wgt_nonadd, ld_width):

    # first snv in 1Mbp
    def _find_first_in_ld_width(wgt, ld_width):
        wgt = wgt.copy()
        wgt.loc[:, 'covered'] = False
        wgt.loc[:, 'first_in_ld'] = False

        while len(wgt.loc[~wgt['covered'], :]) > 0:
            idx_min = wgt.loc[~wgt['covered'], :].index[0]
            print('idx', idx_min)
            wgt.loc[idx_min, 'first_in_ld'] = True

            chrom_min = wgt.loc[idx_min, 'chrom']
            pos_min = wgt.loc[idx_min, 'pos']
            print('chrom,pos', chrom_min, pos_min)
            filt = (wgt['chrom'] == chrom_min) & (wgt['pos'] > pos_min - ld_width) & (wgt['pos'] < pos_min + ld_width)
            wgt.loc[filt, 'covered'] = True
        return wgt

    def _extract_first_in_ld_width(wgt, ld_width):
        wgt = _find_first_in_ld_width(wgt, ld_width)

        wgt = wgt.loc[wgt['first_in_ld'], :]
        wgt = wgt.drop(columns=['covered', 'first_in_ld'])
        return wgt

    # ld_width = 1000000
    wgt_nonadd = _extract_first_in_ld_width(wgt_nonadd, ld_width)
    return wgt_nonadd


def run_gene_catalog(dout, dresult, ddata, mid_path,
                     fpheno, phe, cov_names, fcv, cvi,
                     fwgt, fmap_hg38, nsnv,
                     fgenot, genot_format,
                     ffreq,
                     fcatalog, fgene, fgene_rna):

    method = 'boosting'
    wgt = wgtop.load_wgt_method(fwgt, method)
    wgt = wgtop.format_wgt_boost(wgt, None, '1')
    wgt = format_wgt(wgt)

    logger.debug('wgt {}'.format(wgt))

    use_sid_tmp = True

    # freq
    wgt = merge_freq(wgt, ffreq, use_sid_tmp)

    logger.debug('wgt {}'.format(wgt.head()))
    logger.debug('wgt {}'.format(wgt.info()))

    # TMP: add rs
    if use_sid_tmp:
        wgt = add_rs(wgt, fgenot, genot_format)
        logger.debug('wgt {}'.format(wgt.head()))

    # add gwas
    gmodels = ssop.gmodels_ss()
    ssd = ssop.load_ssd_raw(dresult, phe, cvi, gmodels)
    wgt = ssop.merge_wgt_ssd(wgt, ssd, use_sid_tmp=use_sid_tmp)
    wgt = ssop.judge_genetic_model_gwas(wgt)

    # odds from genot
    # wgt = calculate_odd(wgt, fgenot, genot_format, fpheno, phe, fcv, cov_names)

    # add gene and gene or rna
    wgt = add_gene(wgt, fgene, fgene_rna)
    logger.debug('wgt {}'.format(wgt.info()))

    # add catalog
    wgt = add_catalog(wgt, fcatalog, phe, fmap_hg38, use_sid_tmp)
    # logger.debug('wgt {}'.format(wgt))
    logger.debug('wgt {}'.format(wgt.info()))

    write_paper(dout, wgt, prefix='.all')

    # lastly, extract nsnv
    wgt_nsnv = wgt.loc[wgt['boost_iteration_uniq'] < nsnv, :]
    write_paper(dout, wgt_nsnv, prefix='')

    wgt_nonadd = extract_nonadd_table(wgt_nsnv)
    write_paper(dout, wgt_nonadd, prefix='.nonadd')

    wgt_nonadd_rna = extract_nonadd_table(wgt_nsnv, True)
    write_paper(dout, wgt_nonadd_rna, prefix='.rna.nonadd')

    ld_width = 1000000
    wgt_nonadd_first = extract_wgt_first_in_ld_width(wgt_nonadd, ld_width)
    write_paper(dout, wgt_nonadd_first, prefix='.nonadd.first_in_' + str(ld_width))

    wgt_nonadd_first_rna = extract_wgt_first_in_ld_width(wgt_nonadd_rna, ld_width)
    write_paper(dout, wgt_nonadd_first_rna, prefix='.rna.nonadd.first_in_' + str(ld_width))


def main():
    args = argument()

    run_gene_catalog(pathlib.Path(args.dout),
                     pathlib.Path(args.dresult),
                     pathlib.Path(args.ddata),
                     args.mid_path,
                     pathlib.Path(args.fpheno),
                     args.phe, args.cov_names,
                     pathlib.Path(args.fcv),
                     int(args.cvi),
                     # args.para_cand, args.para_regon, args.para_com,
                     args.fwgt, args.fmap_hg38,
                     int(args.nsnv),
                     args.fgenot, args.genot_format,
                     args.ffreq,
                     args.fcatalog, args.fgene, args.fgene_rna
                     )


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
