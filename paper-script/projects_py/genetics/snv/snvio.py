
from logging import getLogger
import pandas as pd
import numpy as np


logger = getLogger(__name__)


# columns include at least chrom, pos, id, ref, alt
def load_snvs(fsnv, header='infer'):

    def postprocess(snv, has_header):
        if has_header:
            # in case for plink2
            snv.columns = snv.columns.str.replace('#', '')
            col_map = {'CHROM': 'chrom', 'POS': 'pos', 'ID': 'id', 'REF': 'ref', 'ALT': 'alt'}
            snv = snv.rename(columns=col_map)
        else:
            # no header
            if len(snv.columns) == 5:
                snv.columns = 'chrom,pos,id,ref,alt'.split(',')
            else:
                raise RuntimeError('unknown number of columns: ', len(snv.columns), ' in ', fsnv)
        return snv

    if header == 'infer':
        snv = pd.read_csv(fsnv, delim_whitespace=True, dtype=str)
        has_header = (snv.columns.tolist() != np.arange(len(snv.columns)).tolist())
        snv = postprocess(snv, has_header)
    elif header:
        snv = pd.read_csv(fsnv, header=0, delim_whitespace=True, dtype=str)
        snv = postprocess(snv, True)
    else:
        snv = pd.read_csv(fsnv, header=None,
                          delim_whitespace=True, dtype=str)
        snv = postprocess(snv, False)

    return snv


# .bed file created by liftover
def load_hm3(fhm3):
    cols = 'chrom_str,pos,N1,N2'.split(',')
    hm3 = pd.read_csv(fhm3,
                      delim_whitespace=True,
                      header=None, names=cols, dtype=str)

    def f_chrom_str_to_chrom(x):
        return x[3:]

    hm3.loc[:, 'chrom'] = hm3['chrom_str'].map(f_chrom_str_to_chrom)

    hm3 = hm3.drop(columns=['chrom_str'])

    return hm3


def load_fmap_hg38(fmap_hg38):
    names = 'chrom_str,pos,N,id'.split(',')
    df = pd.read_csv(fmap_hg38, delim_whitespace=True, names=names, header=None, dtype=str)
    # dtype=str)

    def fn_chrom_str_to_chrom(x):
        # ex. chr7_KI...
        if '_' not in x:
            return x[3:]
        else:
            return x.split('_')[0][3:]

    df.loc[:, 'chrom'] = df['chrom_str'].map(fn_chrom_str_to_chrom)
    df = df.drop_duplicates(subset=['id'], ignore_index=True)

    df = df.rename(columns={'chrom': 'chrom_hg38', 'pos': 'pos_hg38'})

    # chrom=str
    df = df.astype({'pos_hg38': int})

    return df


#def _prep_gene(gene):
#    # chrom=1-22
#    gene = gene[gene['chrom'].str.isdigit()].copy()
#    gene.loc[:, 'chrom'] = gene['chrom'].map(int)
#
#    gene = gene.sort_values(['chrom', 'start'])
#    gene = gene.reset_index(drop=True)
#
#    print("gene", gene.head(50))
#
#    return gene


def load_gene(fgene,is_rna=False):
    cols = 'chrom,N1,N2,start,end,N3,N4,N5,comment'.split(',')

    gene = pd.read_csv(fgene, sep='\t', header=None, names=cols)

    # chrom=1-22
    gene = gene[gene['chrom'].str.isdigit()].copy()
    # chrom=str
    #gene.loc[:, 'chrom'] = gene['chrom'].map(int)

    gene = gene.sort_values(['chrom', 'start'])
    gene = gene.reset_index(drop=True)

    # logger.debug("gene: {}".format(gene.head(50)))
    logger.debug("gene: {}".format(gene.head()))

    # extract gene name from comment
    def extract_name_from_comment(comment):
        return comment.split(';')[2].strip().split(' ')[1].replace('"', '')

    gene.loc[:, 'gene'] = gene['comment'].map(extract_name_from_comment)
    logger.debug("gene: {}".format(gene[['comment', 'gene']]))

    if is_rna:
        def extract_name_from_comment(comment):
            return comment.split(';')[4].strip().split(' ')[1].replace('"', '')

        gene.loc[:, 'biotype'] = gene['comment'].map(extract_name_from_comment)
        logger.debug("gene biotype {}".format(gene[['comment', 'biotype']]))

    return gene


def load_tsv_header(fin):
    df = pd.read_csv(fin, delim_whitespace=True, header=0, dtype=str)
    return df


def load_assoc(fassoc):
    assoc = load_tsv_header(fassoc)
    return assoc


# TODO:
# def merge_freq(snv, freq):
#    #freq = plink.load_freq_plink2(ffreq)
#    print('freq', freq.info())
#    freq_check = freq[['id', 'ref', 'alt']].set_index('id')
#    snv_check = snv.set_index(('id'))
#    if (freq_check['ref'] != snv_check['ref']).any() or (freq_check['alt'] != snv_check['alt']).any():
#        raise RuntimeError('unmatch allele')
#
#    freq = freq[['id', 'alt_frq']]
#    snv = pd.merge(snv, freq, on='id', how='left')
