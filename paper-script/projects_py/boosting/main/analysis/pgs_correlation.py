
import argparse
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

    parser.add_argument('--phe',
                        help='')

    parser.add_argument('--cvi',
                        help='')

    parser.add_argument('--fwgt',
                        help='')

    parser.add_argument('--fmap_hg38',
                        help='')

    parser.add_argument('--nsnv',
                        help='')

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



def load_catalog(fcatalog,phe):
    phed=phed_catalog()
    phe_catalog=phed[phe]

    catalog = pd.read_csv(fcatalog, sep='\t', header=0, low_memory=False, encoding='utf-8')

    print("catalog", catalog)
    print("catalog", catalog.info())

    catalog = catalog[catalog['DISEASE/TRAIT'].str.lower() == phe_catalog]
    #catalog = catalog[catalog['DISEASE/TRAIT'].str.lower().str.contains(phenotype_catalog)]

    if len(catalog) == 0:
        logger.error('phenotype not found: {}'.format(phe))
        raise RuntimeError('No catalog of disease: ', phe_catalog)

    catalog = catalog.rename(columns={'CHR_ID': 'chrom', 'CHR_POS': 'pos'})

    # exclude exception
    catalog = catalog.dropna(subset=['pos'])
    catalog = catalog[~catalog['pos'].str.contains('x')].copy()

    # pos might be '123;124;125;
    def f_pos(x):
        if ';' not in x:
            return x
        else:
            return x.split(';')[0]
    catalog.loc[:, 'pos'] = catalog['pos'].map(f_pos).astype(int)

    catalog = catalog.astype({'pos': np.float})

    print("catalog", catalog[['FIRST AUTHOR', 'DISEASE/TRAIT', 'chrom', 'pos']])

    if len(catalog)==0:
        raise RuntimeError('No snv in catalog after qc.')


def map_wgt_to_catalog(wgt, catalog):

    wgt.loc[:, 'catalog'] = False
    wgt.loc[:, 'catalog_rel'] = 'NaN'
    wgt.loc[:, 'catalog_gene'] = 'NaN'

    # 1Mbp
    bp_catalog_nearest = 1000000

    for ind in wgt.index:
        # use hg38
        chrom, pos = wgt.loc[ind, ['chrom_hg38', 'pos_hg38']]
        # print("chrom,pos",chrom,pos)

        catalog_chrom = catalog[catalog['chrom'] == str(chrom)]
        # print("catalog",catalog_chrom)

        if len(catalog_chrom) == 0:
            continue

        ind_nearest = catalog_chrom.index[(catalog_chrom['pos'] - pos).abs().argsort()][0]

        diff = pos - int(catalog_chrom.loc[ind_nearest, 'pos'])
        # print("pos",int(catalog_chrom.loc[ind_nearest,'pos']))
        # print("diff",diff)

        if np.abs(diff) < bp_catalog_nearest:
            wgt.loc[ind, 'catalog'] = True
            wgt.loc[ind, 'catalog_rel'] = str(diff) + 'bp' if diff <= 0 else '+' + str(diff) + 'bp'
            wgt.loc[ind, 'catalog_gene'] = str(catalog_chrom.loc[ind_nearest, 'MAPPED_GENE'])
            # wgt.loc[ind,'catalog_gene']=catalog_chrom.loc[ind_nearest,'MAPPED_GENE']

    return wgt


def annotate_gene_overlap(wgt, gene):
    # should allow overlapping genes
    # if snv is in several genes, write all of them
    # if not in any gene, write nearest gene

    bp_gene_nearest = 100000

    wgt.loc[:, 'gene'] = 'NaN'
    wgt.loc[:, 'gene_rel'] = 'NaN'
    if 'biotype' in gene:
        wgt.loc[:, 'gene_biotype'] = 'NaN'

    # regard +- 100kb as 'gene'
    # put 'in' in 'gene_rel' if in the gene
    # put '+123bp'

    for ind in wgt.index:
        chrom, pos = wgt.loc[ind, ['chrom', 'pos']]
        #chrom, pos = wgt.loc[ind, ['chrom', 'pos']]

        gene_chrom = gene[gene['chrom'] == chrom]

        # check if in gene
        gene_in = gene_chrom[(gene_chrom['start'] <= pos) & (pos < gene_chrom['end'])]
        if len(gene_in) > 0:
            #assert len(gene_in) == 1, "snv in several genes"
            wgt.loc[ind, 'gene'] = gene_in['gene'].str.cat(sep=',')
            #wgt.loc[ind, 'gene'] = gene_in.iloc[0, :]['gene']
            wgt.loc[ind, 'gene_rel'] = 'in'
            if 'biotype' in gene:
                wgt.loc[ind, 'gene_biotype'] = gene_in.iloc[0, :]['biotype']
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
                wgt.loc[ind, 'gene'] = gene_chrom.loc[ind_nearest_start, 'gene']
                wgt.loc[ind, 'gene_rel'] = str(diff_start) + 'bp'
                if 'biotype' in gene:
                    wgt.loc[ind, 'gene_biotype'] = gene_chrom.loc[ind_nearest_start, 'biotype']
            else:
                wgt.loc[ind, 'gene'] = gene_chrom.loc[ind_nearest_end, 'gene']
                wgt.loc[ind, 'gene_rel'] = '+' + str(diff_end) + 'bp'
                if 'biotype' in gene:
                    wgt.loc[ind, 'gene_biotype'] = gene_chrom.loc[ind_nearest_end, 'biotype']

    print("wgt", wgt[['sida', 'gene', 'gene_rel']])

    return wgt



def add_catalog(wgt, fcatalog,fmap_hg38,fgene):
    map_hg38=snvio.load_fmap_hg38(fmap_hg38)

    wgt=wgt.merge(map_hg38[['chrom_hg38','pos_hg38','id']],on='id',how='left')

    catalog=load_catalog(fcatalog)

    wgt=map_wgt_to_catalog(wgt,catalog)

    gene=snvio.load_gene(fgene)

    wgt=annotate_gene_overlap(wgt,gene)

    # TODO: fgene_rna


    return wgt


def run_gene_catalog(dout, dresult, ddata, mid_path,
                     phe, cvi,
                     fwgt, fmap_hg38, nsnv,
                     ffreq,
                     fcatalog, fgene, fgene_rna):

    method = 'boosting'
    wgt = wgtop.load_wgt_method(fwgt, method)
    wgt = wgtop.format_wgt_boost(wgt, None, '1')

    logger.debug('wgt {}'.format(wgt))

    use_sid_tmp = True

    # freq
    wgt = merge_freq(wgt, ffreq, use_sid_tmp)

    logger.debug('wgt {}'.format(wgt))
    logger.debug('wgt {}'.format(wgt.info()))

    # TMP: add rs
    if use_sid_tmp:
        wgt = add_rs(wgt, fgenot)

    # lastly, extract nsnv
    # wgt_snv=wgt.iloc[nsnv,:]

    
    


def main():
    args = argument()

    run_gene_catalog(pathlib.Path(args.dout),
                     pathlib.Path(args.dresult),
                     pathlib.Path(args.ddata),
                     args.mid_path,
                     args.phe,
                     int(args.cvi),
                     # args.para_cand, args.para_regon, args.para_com,
                     args.fwgt, args.fmap_hg38,
                     args.nsnv,
                     args.ffreq,
                     args.fcatalog, args.fgene, args.fgene_rna
                     )


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
