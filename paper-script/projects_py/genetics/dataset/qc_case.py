"""
For
- QC on samples or snps.
- split dataset.

input is plink file .fam or .bim.
output is extracted samples file or snps files: "${fout}.samples", "${fout}.snps"

After running this program, you would
`$ plink --bfile $fplink --extract ${fplink}.samples --keep ${fplink}.snps --out ${fplink_out}`


TODO: delete --split_chr and judge by '%' in fgt.
"""

import os
import time

import argparse

import random

import numpy as np
import pandas as pd

import io_gwas as iog
# import validate_ukb


def argument():
    parser = argparse.ArgumentParser(
        prog='gwas',
        usage='gwas',
        description='description',
        epilog='description end',
        add_help=True,
    )

    parser.add_argument('--valid',
                        action='store_true',
                        help="qc on valid samples")

    # with --valid
    parser.add_argument('--only_female',
                        action='store_true',
                        help="only female")

    parser.add_argument('--snps',
                        action='store_true',
                        help="qc on snps file")

    parser.add_argument('--snps_dup_impinfo',
                        action='store_true',
                        help="qc on snps file on duplicated and impinfo")

    parser.add_argument('--snps_hm3',
                        action='store_true',
                        help="")

    parser.add_argument('--snps_hm3_freq',
                        action='store_true',
                        help="")

    parser.add_argument('--down_impinfo',
                        action='store_true',
                        help="down sampling for ldak input")

    parser.add_argument('--down_cal',
                        action='store_true',
                        help="down sampling for ldak input")

    parser.add_argument('--samples_female',
                        action='store_true',
                        help="qc on unrelated samples")

    parser.add_argument('--unrelated',
                        action='store_true',
                        help="qc on unrelated samples")

    parser.add_argument('--unrelated_va',
                        action='store_true',
                        help="qc on unrelated samples")

    parser.add_argument('--unrelated_va12',
                        action='store_true',
                        help="qc on unrelated samples")

    parser.add_argument('--unrelated_test',
                        action='store_true',
                        help="qc on unrelated samples")

    parser.add_argument('--training_propcase',
                        action='store_true',
                        help="training prop case")

    parser.add_argument('--split_covs',
                        action='store_true',
                        help="")

    parser.add_argument('--align_covs',
                        action='store_true',
                        help="")
    parser.add_argument('--align_multiphe',
                        action='store_true',
                        help="")

    parser.add_argument('--zeroone_multiphe',
                        action='store_true',
                        help="")

    parser.add_argument('--cov',
                        action='store_true',
                        help="cov file")

    parser.add_argument('--case',
                        action='store_true',
                        help="case iid")

    parser.add_argument('--case_multi',
                        action='store_true',
                        help="case iid")

    parser.add_argument('--plink_assoc_cov_pcfromukb',
                        action='store_true',
                        help="cov file for plink --assoc")

    parser.add_argument('--plink_assoc_cov_pcfromukb_com2one',
                        action='store_true',
                        help="cov file for plink --assoc")

    parser.add_argument('--plink_assoc_cov',
                        action='store_true',
                        help="cov file for plink --assoc")

    parser.add_argument('--plink_assoc_cov_use_sscore',
                        action='store_true',
                        help="use .sscore.eigenvec instaead of .eigenvec")

    parser.add_argument('--plink_assoc_cov_one2chrom',
                        action='store_true',
                        help="cov file for plink --assoc for each chrom: samples order is different")

    parser.add_argument('--add_a2_plinkassoc',
                        action='store_true',
                        help="add A2 to plinkassoc")

    parser.add_argument('--adjust_pc_projection',
                        action='store_true',
                        help="plink PCA projection does not adjust scale. adjust here.")

    parser.add_argument('--cov_standardize',
                        action='store_true',
                        help='standardize cov for plink2 glm')

    parser.add_argument('--assoc_hm3',
                        action='store_true',
                        help="")

    parser.add_argument('--sida_rs',
                        action='store_true',
                        help="map sida to rs")

    parser.add_argument('--gene_exon',
                        action='store_true',
                        help="")

    parser.add_argument('--snps_hm3_sida_rs',
                        action='store_true',
                        help="map sida to rs")

    parser.add_argument('--fout',
                        default=None,
                        help="output prefix for each chromosome.")

    parser.add_argument('--fout_com',
                        default=None,
                        help="output prefix common for chromosome.")

    parser.add_argument('--fplink',
                        default=None,
                        help='prefix of plink')

    parser.add_argument('--fplink_rs',
                        default=None,
                        help='prefix of plink')

    parser.add_argument('--fphe',
                        default=None,
                        help='phe for multiphe')

    parser.add_argument('--fhm3',
                        default=None,
                        help="")

    parser.add_argument('--fhm3_snvs',
                        default=None,
                        help=".hm3.snvs")

    parser.add_argument('--fsida_rs',
                        default=None,
                        help="")
    parser.add_argument('--fsida_rs_hm3',
                        default=None,
                        help="")

    parser.add_argument('--chrom_split',
                        action='store_true',
                        help='chromosomes are split in fplink')

    parser.add_argument('--phes',
                        type=str,
                        nargs='+',
                        default=None,
                        help='')

    parser.add_argument('--fqc_whitebritish',
                        default=None,
                        help='DataField22006: file of white brithsh')

    parser.add_argument('--fqc_het_missing_outliers',
                        default=None,
                        help='DataField22027: file of het_missing_outliers, samples are to be excluded')

    parser.add_argument('--fqc_sex_chromosome_aneuploidy',
                        default=None,
                        help='DataField22019: file of sex_chromosome_aneuploidy, samples are to be excluded')

    parser.add_argument('--fqc_in_pca',
                        default=None,
                        help='DataField22020: file of in_pca samples')

    parser.add_argument('--fqc_unrelated',
                        default=None,
                        help='DataField22021: file of kinship, samples with 0 (= no kinship was found) are to be used')

    parser.add_argument('--fqc_pc',
                        default=None,
                        help='DataField22009: file of pca')

    parser.add_argument('--fqc_withdrawn',
                        default=None,
                        help='iids withdrawn')

    parser.add_argument('--fsex',
                        default=None,
                        help='DataField31: file of sex')

    parser.add_argument('--fdatebirth',
                        default=None,
                        help='DataField34 + 52: file of year of birth and month of birth')

    parser.add_argument('--fcase',
                        default=None,
                        help='dir with case iid list')

    parser.add_argument('--dcase',
                        default=None,
                        help='for --case_multi')

    parser.add_argument('--fss',
                        default=None,
                        help='gwas summary statistics')

    parser.add_argument('--fss_ori',
                        default=None,
                        help='gwas summary statistics')

    parser.add_argument('--fplinkassoc',
                        default=None,
                        help='plink assoc file')

    parser.add_argument('--fimpinfo',
                        default=None,
                        help='snp impinfo file')

    parser.add_argument('--fcal',
                        default=None,
                        help='snp cal file')

    parser.add_argument('--fmaf',
                        default=None,
                        help='snp maf file to exclude duplicated snps')

    parser.add_argument('--fcov',
                        default=None,
                        help='input cov of sex and age for plink_assoc_cov')

    parser.add_argument('--fcov_out',
                        default=None,
                        help='')

    parser.add_argument('--fphe_out',
                        default=None,
                        help='')

    parser.add_argument('--feigenval',
                        default=None,
                        help='eigenval')

    parser.add_argument('--feigenvec',
                        default=None,
                        help='eigenvec')

    parser.add_argument('--fsscore_training',
                        default=None,
                        help='fsscore_training')

    parser.add_argument('--fsscore',
                        default=None,
                        help='fsscore')

    parser.add_argument('--fgene_exon',
                        default=None,
                        help='fsscore')

    parser.add_argument('--fgene_exon_out',
                        default=None,
                        help='fsscore')

    parser.add_argument('--ffreq',
                        default=None,
                        help='snp maf file to exclude duplicated snps')
    parser.add_argument('--ffreq_hm3',
                        default=None,
                        help='snp maf file to exclude duplicated snps')

    parser.add_argument('--control_ratio',
                        default=1.0,
                        type=float,
                        help='ratio of control/case')

    parser.add_argument('--n_cv',
                        default=5,
                        type=int,
                        help='n_cv')

    parser.add_argument('--seed',
                        type=int,
                        default=None,
                        help='for trainig_prop')

    parser.add_argument('--year_now',
                        default=None,
                        type=int,
                        help='year now for age')

    parser.add_argument('--month_now',
                        default=None,
                        type=int,
                        help='month now for age')

    args = parser.parse_args()

    print(args)

    return args


def read_fam(fplink, chrom_split=False):
    cols = 'fid,iid,N1,N2,sex,phe'.split(',')
    print(cols)
    if chrom_split:
        fplink = fplink.replace('%', '1')
    phe = pd.read_csv(fplink + '.fam', sep='\\s+',
                      header=None,
                      names=cols)
    phe.loc[:, ['fid', 'iid']] = phe[['fid', 'iid']].astype(np.int)
    return phe


def read_fam_str(fplink, chrom_split=False):
    cols = 'fid,iid,N1,N2,sex,phe'.split(',')
    print(cols)
    if chrom_split:
        fplink = fplink.replace('%', '1')
    phe = pd.read_csv(fplink + '.fam', sep='\\s+',
                      dtype=str,
                      header=None,
                      names=cols)
    #phe[['fid', 'iid']] = phe[['fid', 'iid']].astype(np.int)
    return phe


def qc_validiid(fam):
    # exclude invalid negative iid

    fam = fam.loc[fam.index >= 0, :]

    return fam


def qc_include(fam, fqc):
    cols = 'eid,qc'.split(',')
    qc = pd.read_csv(fqc, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')
    qc = qc.sort_index()

    qc = qc[qc['qc'] == 1]

    # here, index name of fam vanished
    index_intersection = fam.index.intersection(qc.index)
    fam = fam.loc[index_intersection]

    return fam


def qc_exclude(fam, fqc):
    cols = 'eid,qc_ex'.split(',')
    qc = pd.read_csv(fqc, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')
    qc = qc.sort_index()

    qc = qc[qc['qc_ex'] != 1]

    index_intersection = fam.index.intersection(qc.index)
    fam = fam.loc[index_intersection]

    return fam


def qc_exclude_iid(fam, fqc):
    cols = 'eid'.split(',')
    qc = pd.read_csv(fqc, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')

    index_intersection = fam.index.intersection(qc.index)

    fam = fam.drop(index=index_intersection)

    return fam


def qc_unrelated(fam, fqc):
    cols = 'eid,qc'.split(',')
    qc = pd.read_csv(fqc, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')
    qc = qc.sort_index()

    # 0 is unrelated
    qc = qc[qc['qc'] == 0]

    # here, index name of fam vanished
    index_intersection = fam.index.intersection(qc.index)
    fam = fam.loc[index_intersection]

    return fam


def qc_female(fam, fsex):
    cols = 'eid,sex'.split(',')
    qc = pd.read_csv(fsex, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')

    fam['sex'] = qc['sex']

    fam = fam[fam['sex'] == 0]

    fam = fam.drop(columns=['sex'])

    return fam


def add_sex(fam, fsex):
    cols = 'eid,sex'.split(',')
    qc = pd.read_csv(fsex, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')

    # fam['sex'] = -9
    # qc = qc[(qc['sex'] == 0) & (qc['sex'] == 1)]
    # index_intersection = fam.index.intersection(qc.index)
    # fam.loc[index_intersection, 'sex'] = qc.loc[index_intersection, 'sex']

    # overwrite 'sex'
    fam['sex'] = qc['sex']

    return fam


def add_age(fam, fdatebirth, year_now, month_now):
    """
    age in 1st, month_now. year_now
    assume all were born in 1st day
    """

    # """
    year_now = 2020
    month_now = 3
    # """

    cols = 'eid,year,month'.split(',')
    qc = pd.read_csv(fdatebirth, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')

    qc['age'] = year_now - qc['year']
    # one more age on those born until this month
    filt = (qc['month'] < month_now)
    qc.loc[filt, 'age'] = qc.loc[filt, 'age'] + 1

    print("age", qc)

    """
    qc['datebirth'] = qc['year'].apply(str) + '-' + qc['month'].apply(str)+'-1'
    qc['now'] = '2020-2-1'

    qc['datebirth'] = pd.to_datetime(
        qc['datebirth'], format='%Y-%m-%d')
    qc['now'] = pd.to_datetime(
        qc['now'], format='%Y-%m-%d')

    print(qc['now'] - qc['datebirth'])

    print("date",qc)

    qc['age'] = (qc['now'] - qc['datebirth']
                 ).apply(lambda v: v.year).astype(int)

    print("qc age", qc)
    """

    fam['age'] = qc['age']

    return fam


def write_index_fam_pd(fout, fam):

    fam[['fid', 'iid']].to_csv(fout, header=None, index=False, sep='\t')


def write_index_fam(fout, fam):
    iid_ = list(fam['iid'].values)
    iid_ = [str(v) + '\t' + str(v) for v in iid_]
    iid_ = '\n'.join(iid_)
    with open(fout, 'w') as f:
        f.write(iid_)


def write_cov(fout, fam):
    fam.to_csv(fout, sep='\t',
               index=False)


def extract_samples_(args):
    fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    fam = fam.sort_index()

    print("fam", fam.head())
    print("bfr all qc", len(fam))
    fam = qc_validiid(fam)
    print("fam", fam.head())
    print("afr qc valid", len(fam))
    fam = qc_exclude_iid(fam, args.fqc_withdrawn)
    print("fam", fam.head())
    print("afr exclude withdrawn", len(fam))
    fam = qc_include(fam, args.fqc_whitebritish)
    print("afr whitebritish", len(fam))
    print("fam", fam.head())
    fam = qc_exclude(fam, args.fqc_sex_chromosome_aneuploidy)
    print("afr sex chromosome aneuploidy", len(fam))
    fam = qc_exclude(fam, args.fqc_het_missing_outliers)
    print("afr het missing outliers", len(fam))

    fam = qc_include(fam, args.fqc_in_pca)
    print("afr in pca", len(fam))

    if args.only_female:
        fam = qc_female(fam, args.fsex)
        print("afr female only", len(fam))

    print("fam", fam.head())
    fam.index.name = 'iid'
    print("fam", fam.head())
    fam = fam.reset_index(drop=False)
    print("fam", fam.head())

    print("output to:", args.fout_com + '.valid.samples')

    write_index_fam(args.fout_com + '.valid.samples', fam)


def split_covs(args):
    # tr.unrelated

    cov = pd.read_csv(args.fcov, sep='\\s+',
                      header=0)
    cov = cov.set_index('iid')
    print('cov', cov)

    for i_cv in range(args.n_cv):
        print("i_cv", i_cv)
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples', sep='\\s+',
                          header=None,
                          names=cols)
        fam.loc[:, ['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
        #fam = fam.set_index('iid')
        #fam = fam.sort_index()
        print("fam", fam.head())

        index_fam = cov.index.intersection(fam['iid'])

        cov_assoc = cov.loc[index_fam]
        cov_assoc = cov_assoc.reset_index()

        col_order = ['fid', 'iid', 'sex', 'age'] + ['PC' + str(i + 1) for i in range(10)]
        cov_assoc = cov_assoc[col_order]

        #print("fam", fam.head())

        #fam.index.name = 'iid'
        #print("fam", fam.head())
        #fam = fam.reset_index(drop=False)
        #print("fam", fam.head())

        #cov_assoc = cov.join(fam, how='inner')
        print("cov_assoc", cov_assoc.head())

        print("output to:", args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.eigenvec_sex_age.cov')

        #write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples', fam)
        write_cov(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.eigenvec_sex_age.cov', cov_assoc)


def align_covs(args):
    # tr.unrelated

    cov = pd.read_csv(args.fcov, sep='\\s+',
                      header=0)
    cov = cov.set_index('iid')
    print('cov', cov)

    cols = 'fid,iid'.split(',')
    fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    # fam = pd.read_csv(args.fsamples, sep='\\s+',
    #                    header=None,
    #                    names=cols)
    fam.loc[:, ['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

    # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    #fam = fam.sort_index()
    print("fam", fam.head())

    index_fam = fam.index.intersection(cov.index)
    #index_fam = cov.index.intersection(fam['iid'])

    cov_assoc = cov.loc[index_fam]
    cov_assoc = cov_assoc.reset_index()

    if 'sex' in cov.columns:
        col_order = ['fid', 'iid', 'sex', 'age'] + ['PC' + str(i + 1) for i in range(10)]
    else:
        col_order = ['fid', 'iid', 'age'] + ['PC' + str(i + 1) for i in range(10)]
    cov_assoc = cov_assoc[col_order]

    #print("fam", fam.head())

    #fam.index.name = 'iid'
    #print("fam", fam.head())
    #fam = fam.reset_index(drop=False)
    #print("fam", fam.head())

    #cov_assoc = cov.join(fam, how='inner')
    print("cov_assoc", cov_assoc.head())

    print("output to:", args.fcov_out)

    #write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples', fam)
    write_cov(args.fcov_out, cov_assoc)


def align_multiphe(args):
    # tr.unrelated

    phe = pd.read_csv(args.fphe, sep='\\s+',
                      dtype=str,
                      header=0)
    cols_phe = phe.columns
    phe = phe.set_index('iid')
    print('phe', phe)
    print('phe', phe.info())

    fam = read_fam_str(args.fplink, chrom_split=args.chrom_split)
    #fam.loc[:, ['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

    # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    #fam = fam.sort_index()
    print("fam", fam)
    print("fam", fam.info())

    index_fam = fam.index
    #index_fam = fam.index.intersection(phe.index)
    #index_fam = cov.index.intersection(fam['iid'])

    cov_assoc = phe.loc[index_fam, :]
    cov_assoc = cov_assoc.reset_index(drop=False)

    cov_assoc = cov_assoc[cols_phe]

    #print("fam", fam.head())

    #fam.index.name = 'iid'
    #print("fam", fam.head())
    #fam = fam.reset_index(drop=False)
    #print("fam", fam.head())

    #cov_assoc = cov.join(fam, how='inner')
    print("cov_assoc", cov_assoc.head())

    print("output to:", args.fcov_out)

    #write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples', fam)
    write_cov(args.fphe_out, cov_assoc)


def zeroone_multiphe(args):
    # tr.unrelated

    phe = pd.read_csv(args.fphe, sep='\\s+',
                      dtype=str,
                      header=0)
    cols_phe = phe.columns
    phe = phe.set_index(['fid', 'iid'])
    print('phe', phe)
    print('phe', phe.info())

    phe = phe.astype(int)
    phe = phe - 1

    cov_assoc = phe

    cov_assoc = cov_assoc.reset_index(drop=False)
    cov_assoc = cov_assoc[cols_phe]

    print("cov_assoc", cov_assoc.head())
    print("output to:", args.fcov_out)

    write_cov(args.fphe_out, cov_assoc)


def extract_samples_female_(args):

    def extract_female(fout, fin, fsex):
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(fin, sep='\\s+',
                          header=None,
                          names=cols)
        fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        fam = fam.set_index('iid')

        print("fam", fam.head())
        print("bfr all qc", len(fam))
        fam = qc_female(fam, fsex)
        print("afr female", len(fam))

        fam.index.name = 'iid'
        fam = fam.reset_index(drop=False)
        print("fam", fam.head())

        print("output to:", fout)
        write_index_fam_pd(fout, fam)

    # TODO: flexible input
    for i_cv in range(args.n_cv):
        print("i_cv", i_cv)
        extract_female(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.female.samples',
                       args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples',
                       args.fsex)

        extract_female(args.fout_com + '.obs.cv' + str(i_cv) + '.va.unrelated.female.samples',
                       args.fout_com + '.obs.cv' + str(i_cv) + '.va.unrelated.samples',
                       args.fsex)

    # obs,test

    extract_female(args.fout_com + '.obs.unrelated.female.samples',
                   args.fout_com + '.obs.unrelated.samples',
                   args.fsex)

    extract_female(args.fout_com + '.test.unrelated.female.samples',
                   args.fout_com + '.test.unrelated.samples',
                   args.fsex)


def extract_samples_unrelated_(args):
    # for cv0

    # TODO: flexible input
    for i_cv in range(args.n_cv):
        print("i_cv", i_cv)
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.samples', sep='\\s+',
                          header=None,
                          names=cols)
        fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
        fam = fam.set_index('iid')
        fam = fam.sort_index()

        print("fam", fam.head())
        print("bfr all qc", len(fam))
        fam = qc_unrelated(fam, args.fqc_unrelated)
        print("afr unrelated", len(fam))

        print("fam", fam.head())
        fam.index.name = 'iid'
        print("fam", fam.head())
        fam = fam.reset_index(drop=False)
        print("fam", fam.head())

        print("output to:", args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples')

        write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples', fam)

    # for .obs
    cols = 'fid,iid'.split(',')
    fam = pd.read_csv(args.fout_com + '.obs.samples', sep='\\s+',
                      header=None,
                      names=cols)
    fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

    # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    fam = fam.sort_index()

    print("fam", fam.head())
    print("bfr all qc", len(fam))
    fam = qc_unrelated(fam, args.fqc_unrelated)
    print("afr unrelated", len(fam))

    print("fam", fam.head())
    fam.index.name = 'iid'
    print("fam", fam.head())
    fam = fam.reset_index(drop=False)
    print("fam", fam.head())

    print("output to:", args.fout_com + '.obs.unrelated.samples')

    write_index_fam(args.fout_com + '.obs.unrelated.samples', fam)


def extract_samples_unrelated_va_(args):

    # TODO: flexible input
    for i_cv in range(args.n_cv):
        print("i_cv", i_cv)
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(args.fout_com + '.obs.cv' + str(i_cv) + '.va.samples', sep='\\s+',
                          header=None,
                          names=cols)
        fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
        fam = fam.set_index('iid')
        fam = fam.sort_index()

        print("fam", fam.head())
        print("bfr all qc", len(fam))
        fam = qc_unrelated(fam, args.fqc_unrelated)
        print("afr unrelated", len(fam))

        print("fam", fam.head())
        fam.index.name = 'iid'
        print("fam", fam.head())
        fam = fam.reset_index(drop=False)
        print("fam", fam.head())

        print("output to:", args.fout_com + '.obs.cv' + str(i_cv) + '.va.unrelated.samples')

        write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.va.unrelated.samples', fam)


def extract_samples_unrelated_va12_(args):

    # TODO: flexible input
    for i_cv in range(args.n_cv):
        print("i_cv", i_cv)
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(args.fout_com + '.obs.cv' + str(i_cv) + '.va1.samples', sep='\\s+',
                          header=None,
                          names=cols)
        fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
        fam = fam.set_index('iid')
        fam = fam.sort_index()

        print("fam", fam.head())
        print("bfr all qc", len(fam))
        fam = qc_unrelated(fam, args.fqc_unrelated)
        print("afr unrelated", len(fam))

        print("fam", fam.head())
        fam.index.name = 'iid'
        print("fam", fam.head())
        fam = fam.reset_index(drop=False)
        print("fam", fam.head())

        print("output to:", args.fout_com + '.obs.cv' + str(i_cv) + '.va1.unrelated.samples')

        write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.va1.unrelated.samples', fam)

    # TODO: flexible input
    for i_cv in range(args.n_cv):
        print("i_cv", i_cv)
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(args.fout_com + '.obs.cv' + str(i_cv) + '.va2.samples', sep='\\s+',
                          header=None,
                          names=cols)
        fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
        fam = fam.set_index('iid')
        fam = fam.sort_index()

        print("fam", fam.head())
        print("bfr all qc", len(fam))
        fam = qc_unrelated(fam, args.fqc_unrelated)
        print("afr unrelated", len(fam))

        print("fam", fam.head())
        fam.index.name = 'iid'
        print("fam", fam.head())
        fam = fam.reset_index(drop=False)
        print("fam", fam.head())

        print("output to:", args.fout_com + '.obs.cv' + str(i_cv) + '.va2.unrelated.samples')

        write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.va2.unrelated.samples', fam)


def extract_samples_unrelated_test_(args):

    # TODO: flexible input
    cols = 'fid,iid'.split(',')
    fam = pd.read_csv(args.fout_com + '.test.samples', sep='\\s+',
                      header=None,
                      names=cols)
    fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

    # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    fam = fam.sort_index()

    print("fam", fam.head())
    print("bfr all qc", len(fam))
    fam = qc_unrelated(fam, args.fqc_unrelated)
    print("afr unrelated", len(fam))

    print("fam", fam.head())
    fam.index.name = 'iid'
    print("fam", fam.head())
    fam = fam.reset_index(drop=False)
    print("fam", fam.head())

    print("output to:", args.fout_com + '.test.unrelated.samples')

    write_index_fam(args.fout_com + '.test.unrelated.samples', fam)


def extract_controlratio(fam, control_ratio=1.0, seed=None):
    """
    cols = 'eid,qc'.split(',')
    qc = pd.read_csv(fqc, sep=",",
                     names=cols,
                     header=0)
    qc = qc.set_index('eid')
    qc = qc.sort_index()

    # 0 is unrelated
    qc = qc[qc['qc'] == 0]

    # here, index name of fam vanished
    index_intersection = fam.index.intersection(qc.index)
    fam = fam.loc[index_intersection]
    """
    pass

    ncase = len(fam[fam['phe'] == 1])
    ncont = int(ncase * control_ratio)

    print("ncase, ncont", ncase, ncont)

    ncont_in = len(fam[fam['phe'] == 0])

    print("ncont_in", ncont_in)

    np.random.seed(seed)
    cont_iindexs = np.random.choice(ncont_in, ncont, replace=False)

    print("cont_iindexs", cont_iindexs[:5])

    indexs_cont = fam[fam['phe'] == 0].iloc[cont_iindexs].index
    indexs_case = fam[fam['phe'] == 1].index

    indexs = np.append(indexs_cont, indexs_case)

    print("indexs", len(indexs))

    fam = fam.loc[indexs]

    return fam


def extract_training_propcase_(args):
    """
    subsampling control in training dataset in proportion to args.prop_case
    """
    # for cv0

    print("control ratio", args.control_ratio)

    for i_cv in range(args.n_cv):
        # TODO: flexible input
        cols = 'fid,iid'.split(',')
        fam = pd.read_csv(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.samples', sep='\\s+',
                          header=None,
                          names=cols)
        fam[['fid', 'iid']] = fam[['fid', 'iid']].astype(np.int)

        # fam = read_fam(args.fplink, chrom_split=args.chrom_split)
        fam = fam.set_index('iid')
        fam = fam.sort_index()

        # add phe
        cols = 'fid,iid,N1,N2,sex,phe'.split(',')
        fam_phe = read_fam(args.fplink, chrom_split=args.chrom_split)
        fam_phe[['fid', 'iid']] = fam_phe[['fid', 'iid']].astype(np.int)
        fam_phe = fam_phe.set_index('iid')
        fam_phe['phe'] = fam_phe['phe'].astype(np.int)
        fam_phe['phe'] = fam_phe['phe'] - 1

        fam['phe'] = fam_phe.loc[fam.index, 'phe']

        print("fam with phe", fam.head())

        print("fam", fam.head())
        print("bfr all qc", len(fam))
        fam = extract_controlratio(fam, args.control_ratio, args.seed)
        print("afr unrelated", len(fam))
        fam = fam.sort_index()

        print("fam", fam.head())
        fam.index.name = 'iid'
        print("fam", fam.head())
        fam = fam.reset_index(drop=False)
        print("fam", fam.head())

        print("output to:", args.fout_com + '.obs.cv' + str(i_cv) + '.tr.prop.samples')

        write_index_fam(args.fout_com + '.obs.cv' + str(i_cv) + '.tr.prop.samples', fam)


def write_cov_(args):
    fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.drop(columns='N1,N2,sex,phe'.split(','))
    fam = fam.set_index('iid')
    #fam = fam.sort_index()

    fam = add_sex(fam, args.fsex)
    print("fam sex 1=", len(fam[fam['sex'] == 1]),
          ", 0=", len(fam[fam['sex'] == 0]))
    fam = add_age(fam, args.fdatebirth, args.year_now, args.month_now)
    print("fam age max=", fam['age'].max(),
          ", min=", fam['age'].min())

    fam = fam.reset_index(drop=False)
    fam['fid'] = fam['iid']
    cols = 'fid,iid,sex,age'.split(',')
    fam = fam[cols]
    fam[['sex', 'age']] = fam[['sex', 'age']].fillna(-1)
    fam = fam.astype(int)

    print("output to:", args.fout_com + '.cov')

    write_cov(args.fout_com + '.cov', fam)


def set_case(fam, dcase):
    fam['phe'] = 0
    for d in os.listdir(dcase):
        print("dir file", d)
        for f in os.listdir(os.path.join(dcase, d)):
            fname = os.path.join(dcase, d, f)
            if not os.path.isfile(fname):
                continue

            print("fname of case", fname)
            cs = pd.read_csv(fname, sep='\t', header=0)
            cs = cs.set_index('eid')
            cs.index = cs.index.astype(int)
            print("csindex", len(cs.index), "unique", len(cs.index.unique()))

            # when some of phe are not in cs
            phe_case_index = cs.index.unique().intersection(fam.index)

            print("pheindex", len(fam.index),
                  "intersection", len(phe_case_index))

            fam.loc[phe_case_index, 'phe'] = 1
    return fam


def write_case_(args):
    fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    fam = set_case(fam, args.fcase)
    print("case: ", len(fam[fam['phe'] == 1]), "/", len(fam))
    fam = fam[fam['phe'] == 1]
    print("output to:", args.fout_com + '.case')
    fam = fam.reset_index()
    write_index_fam(args.fout_com + '.case', fam)
    # overwrite_fam(args.fplink, fam)


def load_case_multi(file):

    case = pd.read_csv(file, sep='\t',
                       header=0,
                       dtype=str,
                       names=['eid'])
    return case


def write_case_multi_(args):
    fam = read_fam_str(args.fplink, chrom_split=args.chrom_split)
    fam = fam.set_index('iid')
    phes = args.phes
    for phe in phes:
        fcase = args.dcase + phe + '/case.txt'
        case = load_case_multi(fcase)
        print('case', case)
        common_case = case.set_index('eid').index.intersection(fam.index)
        # 1 is control
        fam.loc[:, phe] = '1'
        fam.loc[common_case, phe] = '2'
        # print('fam',fam)
        print("case", phe, ": ", len(fam[fam[phe] == '2']), "/", len(fam))

    print("output to:", args.fout_com + '.phe')
    fam = fam.reset_index(drop=False)
    cols = 'fid,iid'.split(',') + phes
    fam = fam[cols]
    fam.to_csv(args.fout_com + '.phe', sep='\t', index=False)

    # overwrite_fam(args.fplink, fam)


def read_bim(fplink, chrom_split=False, chrom=None):
    cols = 'chrom,vid,N1,pos,A1,A2'.split(',')
    if chrom_split:
        fplink = fplink.replace('%', str(chrom))
    bim = pd.read_csv(fplink + '.bim', sep='\\s+',
                      header=None,
                      names=cols)
    return bim


def extract_snps_(args):
    ss = iog.load_ss(args.fss, nss=None, ssformat='CUSTOM',
                     ss_columns=None, loadonly=True)
    print("ss", len(ss))

    nsnps_all = 0
    for chrom in range(1, 23):
        print("chrom", chrom)
        bim = read_bim(args.fplink, chrom_split=True, chrom=chrom)
        print("original bim", len(bim))

        bim = iog.add_sids(bim)
        bim = iog.sort_genot_sida(bim)

        filter_bim = iog.isin_genot_sida(bim, ss, both=False)
        bim = bim.loc[filter_bim]

        print("isin bim", len(bim))

        snps_inter = list(bim['vid'])

        nsnps_all += len(snps_inter)

        snps_inter_ = '\n'.join(snps_inter)

        fout = args.fout.replace('%', str(chrom))
        print("output to:", fout + '.ss.snps')
        with open(fout + '.ss.snps', 'w') as f:
            f.write(snps_inter_)

    print("number of extracted snps: ", nsnps_all)


def extract_snps_dup_impinfo_(args):
    """
    ss = iog.load_ss(args.fss, nss=None, ssformat='CUSTOM',
                     ss_columns=None, loadonly=True)
    print("ss", len(ss))
    """

    nsnps_all = 0
    nsnps_nondup_all = 0
    nsnps_nondup_info3_all = 0

    for chrom in range(1, 23):
        # exclude duplicated snvs, leave larger maf
        # assume bim and maf are aligned
        # load bim
        print("chrom", chrom)
        bim = read_bim(args.fplink, chrom_split=True, chrom=chrom)
        print("original bim", len(bim))
        bim = iog.add_sids(bim)
        # do not sort here to align frq

        nsnps_all += len(bim)

        # no pos included
        print("load maf")
        cols = 'chrom,vid,A1,A2,maf,N1'.split(',')
        fmaf = args.fmaf.replace('%', str(chrom))
        maf = pd.read_csv(fmaf, sep='\\s+',
                          header=0,
                          names=cols)
        print("maf", maf.head())

        filt_nondup = ~bim['sid'].duplicated(keep=False)

        bim_dup_filt = bim['sid'].duplicated(keep='first')
        bim_dup = bim[bim_dup_filt]
        print("bim sid dup", len(bim_dup))

        print("bim_dup", bim_dup.head())

        for bim_dup_index in bim_dup.index:
            sid = bim_dup.loc[bim_dup_index, 'sid']
            indexs_dup = bim[bim['sid'] == sid].index
            if len(indexs_dup) > 3:
                print("dup seems too many:", sid, indexs_dup)

            print("dup bim", bim.loc[indexs_dup])

            index_max = maf.loc[indexs_dup, 'maf'].idxmax()

            print("index_max", index_max)

            filt_nondup[index_max] = True

        bim = bim[filt_nondup]
        bim = bim.reset_index()
        print("bim nondup", len(bim))

        nsnps_nondup_all += len(bim)

        bim = iog.sort_genot_sida(bim)

        # load impinfo
        # N1 = maf
        print("load impinfo")
        cols = 'vid,rs,pos,A1,A2,N1,N2,info'.split(',')
        finfo = args.fimpinfo.replace('%', str(chrom))
        info = pd.read_csv(finfo, sep='\\s+',
                           header=None,
                           names=cols)
        info['chrom'] = chrom
        info = iog.add_sids(info)
        info = iog.sort_genot_sida(info)

        print("all info", len(info))
        info = info[info['info'] > 0.3]
        print("info>0.3", len(info))

        filter_bim = iog.isin_genot_sida(bim, info, both=False)
        bim = bim.loc[filter_bim]

        print("isin bim", len(bim))

        nsnps_nondup_info3_all += len(bim)

        snps_inter = list(bim['vid'])

        snps_inter_ = '\n'.join(snps_inter)

        fout = args.fout.replace('%', str(chrom))
        print("output to:", fout + '.qcinfo.snps')
        with open(fout + '.qcinfo.snps', 'w') as f:
            f.write(snps_inter_)

    print("number of init snps: ", nsnps_all)
    print("number of nondup snps: ", nsnps_nondup_all)
    print("number of nondup and info>0.3 snps: ", nsnps_nondup_info3_all)


def extract_snps_hm3(args):
    """
    output list of vid of intersection of hm3 and fgt
    """

    fhm3 = args.fhm3

    #cols = 'chrom,vid,N1,pos'.split(',')
    # hm3 = pd.read_csv(fhm3, sep='\\s+',
    #                       header=None,
    #                       names=cols)
    # hm3=hm3.set_index(['chrom','pos'])

    # use liftover
    cols = 'chrom_str,pos,N,vid'.split(',')
    hm3 = pd.read_csv(fhm3, sep='\t', header=None, names=cols)

    def f_chrom_str_to_chrom(x):
        return x[3:]

    hm3.loc[:, 'chrom'] = hm3['chrom_str'].map(f_chrom_str_to_chrom).astype(int)

    hm3 = hm3.set_index(['chrom', 'pos'])

    hm3.loc[:, 'use'] = False
    hm3.loc[:, 'sida'] = ''

    for chrom in range(1, 23):
        # exclude duplicated snvs, leave larger maf
        # assume bim and maf are aligned
        # load bim
        print("chrom", chrom)
        bim = read_bim(args.fplink, chrom_split=True, chrom=chrom)
        print("original bim", len(bim))
        bim = iog.add_sids(bim)
        # do not sort here to align frq

        bim = bim.set_index(['chrom', 'pos'])

        inter_index = bim.index.intersection(hm3.index)

        print("intersect", len(inter_index))

        hm3.loc[inter_index, 'use'] = True
        hm3.loc[inter_index, 'sida'] = bim.loc[inter_index, 'sida']

    hm3 = hm3.loc[hm3['use'], 'sida']
    hm3 = hm3.drop_duplicates()

    hm3.to_csv(args.fout, index=False, header=False)


def extract_assoc_hm3(args):
    """
    output list of vid of intersection of hm3 and fgt
    """

    fhm3 = args.fhm3_snvs

    cols = 'sida'.split(',')
    hm3 = pd.read_csv(fhm3, sep='\t', header=None, names=cols)

    # loading takes time

    ss_ori = pd.read_csv(args.fss_ori, sep='\t', dtype=str)
    # ss_ori=pd.read_csv(args.fss_ori,sep='\t',dtype={'p':str})
    # now p=str is unnecessary since make p<1e-300 in .ss.format
    # -> should not do this
    # ss_ori=pd.read_csv(args.fss_ori,sep='\t',dtype={'p':str})
    # p could be 1e-341 in psr
    # no effect; float64 is at most what pandas can use
    # ss_ori=pd.read_csv(args.fss_ori,sep='\t',dtype={'p':np.longdouble})
    # ss_ori=pd.read_csv(args.fss_ori,sep='\t',dtype={'p':np.float128})
    print('loaded hm3, ss_ori')

    # print(ss_ori['p'].info())
    #print('ss_ori small p',ss_ori.loc[ss_ori['p']<10**(-300)])

    # time almost same? bet. complicated vs simple

    # complicated
    # hm3=hm3.set_index('sida')
    # ss_ori=ss_ori.set_index('rs')
    # index_common=ss_ori.index.intersection(hm3.index)
    # ss=ss_ori.loc[index_common,:]
    # ss=ss.reset_index(drop=False)

    # simple
    ss = pd.merge(ss_ori, hm3[['sida']], left_on='rs', right_on='sida', how='right')

    # drop 'sida'
    ss = ss.drop(columns=['sida'])

    assert len(ss.columns) == len(ss.columns)

    ss.to_csv(args.fss, sep='\t', index=False)


def map_sida_rs(args):

    bims = []
    for chrom in range(1, 23):
        # exclude duplicated snvs, leave larger maf
        # assume bim and maf are aligned
        # load bim
        print("chrom", chrom)
        bim = read_bim(args.fplink, chrom_split=True, chrom=chrom)
        print("original bim", len(bim))
        #bim = iog.add_sids(bim)
        # do not sort here to align frq

        bim = bim.set_index(['chrom', 'pos'])

        bim_rs = read_bim(args.fplink_rs, chrom_split=True, chrom=chrom)
        print("original bim_rs", len(bim_rs))
        #bim = iog.add_sids(bim)

        bim_rs = bim_rs.set_index(['chrom', 'pos'])

        inter_index = bim.index.intersection(bim_rs.index)

        # assume bim is in bim_rs
        assert len(bim) == len(bim.loc[inter_index, :])
        bim_rs = bim_rs.loc[inter_index, :]
        print("bim_rs intersect", len(bim_rs))

        bim = bim.reset_index(drop=False)
        bim = iog.add_sids(bim)

        bim_rs = bim_rs.reset_index(drop=False)
        #bim_rs = iog.add_sids(bim_rs)
        bim_rs.loc[:, 'sida'] = bim_rs['chrom'].astype(str) + ':' + bim_rs['pos'].astype(str) + ':' + bim_rs['A1'] + ':' + bim_rs['A2']
        bim_rs.loc[:, 'sida2'] = bim_rs['chrom'].astype(str) + ':' + bim_rs['pos'].astype(str) + ':' + bim_rs['A2'] + ':' + bim_rs['A1']

        filt = bim_rs['sida'].isin(bim['sida'])
        filt2 = bim_rs['sida2'].isin(bim['sida'])
        filt_both = filt | filt2
        bim_rs = bim_rs.loc[filt_both, :]

        bim = bim.set_index(['chrom', 'pos'])
        bim_rs = bim_rs.set_index(['chrom', 'pos'])
        print('dup', bim_rs.loc[bim_rs.index.duplicated(), :])
        assert len(bim_rs.loc[bim_rs.index.duplicated(), :]) == 0, 'duplicated sid'

        assert len(bim) == len(bim_rs), 'same'
        # bim.loc[inter_index,'sida']=bim_rs.loc[inter_index,'vid']
        bim.loc[:, 'rs'] = bim_rs.loc[:, 'vid']

        bim = bim.reset_index(drop=False)
        cols = 'chrom,pos,A1,A2,vid,rs'.split(',')
        #print('output to',args.fplink.replace('%', str(chrom))+'.sida_rs.tsv')
        #bim[cols].to_csv(args.fplink.replace('%', str(chrom))+'.sida_rs.tsv', sep='\t', index=False)

        bims.append(bim)

    bim = pd.concat(bims, ignore_index=True)
    print('output to', args.fout_com + '.sida_rs.tsv')
    bim.to_csv(args.fout_com + '.sida_rs.tsv', sep='\t', index=False)


def extract_snps_hm3_sida_rs(args):
    """
    output list of vid of intersection of hm3 and fgt
    """

    fhm3 = args.fhm3_snvs

    cols = 'sida'.split(',')
    hm3 = pd.read_csv(fhm3, sep='\t', header=None, names=cols)

    # loading takes time

    freq = pd.read_csv(args.fsida_rs, delim_whitespace=True)
    # ss_ori=pd.read_csv(args.fss_ori,sep='\t')
    print('loaded hm3, ss_ori')

    # print(ss_ori['p'].info())
    #print('ss_ori small p',ss_ori.loc[ss_ori['p']<10**(-300)])

    # time almost same? bet. complicated vs simple

    # complicated
    # hm3=hm3.set_index('sida')
    # ss_ori=ss_ori.set_index('rs')
    # index_common=ss_ori.index.intersection(hm3.index)
    # ss=ss_ori.loc[index_common,:]
    # ss=ss.reset_index(drop=False)

    # simple
    freq_hm3 = pd.merge(freq, hm3[['sida']], left_on='sida', right_on='sida', how='right')

    # drop 'sida'
    # freq_hm3=freq_hm3.drop(columns=['sida'])

    assert len(freq_hm3.columns) == len(freq.columns)

    freq_hm3.to_csv(args.fsida_rs_hm3, sep='\t', index=False)


def extract_snps_hm3_freq(args):
    """
    output list of vid of intersection of hm3 and fgt
    """

    fhm3 = args.fhm3_snvs

    cols = 'sida'.split(',')
    hm3 = pd.read_csv(fhm3, sep='\t', header=None, names=cols)

    # loading takes time

    freq = pd.read_csv(args.ffreq, delim_whitespace=True)
    # ss_ori=pd.read_csv(args.fss_ori,sep='\t')
    print('loaded hm3, ss_ori')

    # print(ss_ori['p'].info())
    #print('ss_ori small p',ss_ori.loc[ss_ori['p']<10**(-300)])

    # time almost same? bet. complicated vs simple

    # complicated
    # hm3=hm3.set_index('sida')
    # ss_ori=ss_ori.set_index('rs')
    # index_common=ss_ori.index.intersection(hm3.index)
    # ss=ss_ori.loc[index_common,:]
    # ss=ss.reset_index(drop=False)

    # simple
    freq_hm3 = pd.merge(freq, hm3[['sida']], left_on='SNP', right_on='sida', how='right')

    # drop 'sida'
    freq_hm3 = freq_hm3.drop(columns=['sida'])

    assert len(freq_hm3.columns) == len(freq.columns)

    freq_hm3.to_csv(args.ffreq_hm3, sep='\t', index=False)


def extract_snps_down_impinfo(args):
    # do not use
    raise NotImplementedError()

    # exclude duplicated snvs, leave larger maf
    # assume bim and maf are aligned
    # load bim
    bim = read_bim(args.fplink, chrom_split=False)
    print("original bim", len(bim))
    bim = iog.add_sids(bim)
    # do not sort here to align frq

    bim = iog.sort_genot_sida(bim)
    print("number of init snps: ", len(bim))

    snps_inter = []

    for chrom in range(1, 23):
        print("chrom", chrom)

        # load impinfo
        # N1 = maf
        print("load impinfo")
        cols = 'vid,rs,pos,A1,A2,N1,N2,info'.split(',')
        finfo = args.fimpinfo.replace('%', str(chrom))
        info = pd.read_csv(finfo, sep='\\s+',
                           header=None,
                           names=cols)
        info['chrom'] = chrom
        info = iog.add_sids(info)
        info = iog.sort_genot_sida(info)

        print("all info", len(info))
        info = info[info['info'] > 0.9]
        print("info>0.9", len(info))

        filter_bim = iog.isin_genot_sida(bim, info, both=False)
        bim_chrom = bim.loc[filter_bim]

        print("isin bim", len(bim_chrom))

        snps_inter_ = list(bim_chrom['vid'])
        snps_inter += snps_inter_

    print("number of info>0.9 snps: ", len(snps_inter))

    snps_inter_ = '\n'.join(snps_inter)

    fout = args.fout
    print("output to:", fout + '.small.qcinfo.snps')
    with open(fout + '.small.qcinfo.snps', 'w') as f:
        f.write(snps_inter_)


def prep_gene(gene):
    # chrom=1-22
    gene = gene[gene['chrom'].str.isdigit()].copy()
    gene.loc[:, 'chrom'] = gene['chrom'].map(int)

    gene = gene.sort_values(['chrom', 'start'])
    gene = gene.reset_index(drop=True)

    print("gene", gene.head(50))

    return gene


def extract_gene_name(gene):
    # extract gene_name and gene_id
    def extract_name_from_comment(comment):
        return comment.split(';')[2].strip().split(' ')[1].replace('"', '')

    def extract_id_from_comment(comment):
        return comment.split(';')[0].strip().split(' ')[1].replace('"', '')

    filt_gene = gene['gene_type'] == 'gene'

    gene.loc[filt_gene, 'gene'] = gene.loc[filt_gene, 'comment'].map(extract_name_from_comment)
    gene.loc[filt_gene, 'gene_id'] = gene.loc[filt_gene, 'comment'].map(extract_id_from_comment)
    print("gene", gene[['comment', 'gene']])

    # this is non-zero so use gene_id
    print('Duplicated gene_name', gene.loc[filt_gene & gene.loc[filt_gene, 'gene'].duplicated(keep=False)].to_string())

    return gene


def extract_gene_name_exon(gene):
    def extract_name_from_comment(comment):
        return comment.split(';')[5].strip().split(' ')[1].replace('"', '')

    def extract_id_from_comment(comment):
        return comment.split(';')[0].strip().split(' ')[1].replace('"', '')

    filt_gene = gene['gene_type'] == 'exon'

    gene.loc[filt_gene, 'gene'] = gene.loc[filt_gene, 'comment'].map(extract_name_from_comment)
    gene.loc[filt_gene, 'gene_id'] = gene.loc[filt_gene, 'comment'].map(extract_id_from_comment)
    print("gene", gene[['comment', 'gene']])

    return gene


def create_gene_exon(args):
    """
    *** USE GENE_ID not GENE_NAME since gene_name are not unique  ***
    """
    fgene_exon = args.fgene_exon
    fout = args.fgene_exon_out

    # ref boosting_cpp_observe_wgt_find_interesting_snv.py

    cols = 'chrom,N1,gene_type,start,end,N3,N4,N5,comment'.split(',')
    usecols = 'chrom,gene_type,start,end,comment'.split(',')
    gene_exon = pd.read_csv(fgene_exon, sep='\t', header=None, names=cols,
                            usecols=usecols, dtype={'chrom': str})
    # row: gene_whole or exon
    # all genes in exon should have gene_whole row and the region should be included by gene_whole

    gene_exon = prep_gene(gene_exon)

    gene_exon.loc[:, 'gene'] = ''
    gene_exon.loc[:, 'gene_id'] = ''
    gene_exon = extract_gene_name(gene_exon)
    gene_exon = extract_gene_name_exon(gene_exon)

    assert len(gene_exon[gene_exon['gene'] == '']) == 0

    gene_exon = gene_exon.drop(columns=['comment'])

    # check all exon region are in gene
    for gene_id in gene_exon['gene_id'].unique():
        # print('gene_name',gene_name)

        gene_ = gene_exon[gene_exon['gene_id'] == gene_id]

        # check gene_whole in unique
        if len(gene_[gene_['gene_type'] == 'gene']) == 0:
            print('Exon region exists but gene row does not exist so exclude this gene.', gene_id)
            print('Exon region exists but gene row does not exist.', gene_)
            # exclude this gene
            gene_exon = gene_exon[gene_exon['gene_id'] != gene_id].copy()
            continue

        if len(gene_[gene_['gene_type'] == 'gene']) >= 2:
            print('Gene rows are duplicated.', gene_)
            raise RuntimeError()

        gene_start = gene_.loc[gene_['gene_type'] == 'gene', 'start'].iloc[0]
        gene_end = gene_.loc[gene_['gene_type'] == 'gene', 'end'].iloc[0]
        for exon_index, exon_row in gene_[gene_['gene_type'] == 'exon'].iterrows():
            if exon_row['start'] < gene_start or exon_row['end'] > gene_end:
                print('Exon region are out of gene region', gene_)
                raise RuntimeError()

    print('No promblem!!')

    gene_exon.to_csv(fout, sep='\t', index=False)


def extract_snps_down_cal(args):

    # exclude duplicated snvs, leave larger maf
    # assume bim and maf are aligned
    # load bim
    bim = read_bim(args.fplink, chrom_split=False)
    print("original bim", len(bim))
    bim = iog.add_sids(bim)
    # do not sort here to align frq

    bim = iog.sort_genot_sida(bim)
    print("number of init snps: ", len(bim))

    snps_inter = []

    for chrom in range(1, 23):
        print("chrom", chrom)

        # load cal
        print("load cal")
        #cols = 'vid,rs,pos,A1,A2,N1,N2,info'.split(',')
        cols = 'chrom,rs,N1,pos,A1,A2'.split(',')
        fcal = args.fcal.replace('%', str(chrom)) + '.bim'
        cal = pd.read_csv(fcal, sep='\\s+',
                          header=None,
                          names=cols)
        cal = iog.add_sids(cal)
        cal = iog.sort_genot_sida(cal)

        filter_bim = iog.isin_genot_sida(bim, cal, both=False)
        #print("filt", np.count_nonzero(filter_bim))

        bim_chrom = bim.loc[filter_bim]
        print("isin bim", len(bim_chrom))

        """
        # no use singe isin_genot_sida contains this
        #snps_inter_ = list(bim_chrom['vid'])
        #snps_inter+=snps_inter_

        # also add flipped ver.
        def add_sida_flip(df):
            df['sida']=df['chrom'].apply(str)+':'+df['pos'].apply(str)+':'+df['A2']+':'+df['A1']
            return df

        cal = add_sida_flip(cal)
        cal = iog.sort_genot_sida(cal)

        filter_bim2 = iog.isin_genot_sida(bim, cal, both=False)
        print("filt2", np.count_nonzero(filter_bim2))

        filt=filter_bim|filter_bim2

        bim_chrom = bim.loc[filt]

        print("isin bim", len(bim_chrom))
        """

        snps_inter_ = list(bim_chrom['vid'])
        snps_inter += snps_inter_

    print("number of info>0.9 snps: ", len(snps_inter))

    snps_inter_ = '\n'.join(snps_inter)

    fout = args.fout
    print("output to:", fout + '.small.cal.snps')
    with open(fout + '.small.cal.snps', 'w') as f:
        f.write(snps_inter_)


def read_eigenvec(fplink, fuse_sscore=False):
    # cols = 'chrom,vid,N1,pos,A1,A2'.split(',')
    if not fuse_sscore:
        eigenvec = pd.read_csv(fplink + '.eigenvec', sep='\\s+',
                               header=0)
    else:
        eigenvec = pd.read_csv(fplink + '.sscore.eigenvec', sep='\\s+',
                               header=0)

    print("eigenvec", eigenvec.head())
    # 'fid' could be read as '#fid'
    eigenvec = eigenvec.rename(columns={'#FID': 'fid', 'IID': 'iid'})

    eigenvec[['fid', 'iid']] = eigenvec[['fid', 'iid']].astype(np.int)
    return eigenvec


def load_qc_ukb(fpc):

    pc_max = 10

    #cols = 'eid,qc'.split(',')
    qc = pd.read_csv(fpc, sep=",",
                     # names=cols,
                     header=0)
    print('qc', qc)
    qc = qc.rename(columns={'eid': 'iid'})
    qc = qc.set_index('iid')
    print('qc', qc)
    qc.columns = ['PC' + str(i + 1) for i, col in enumerate(qc.columns)]
    qc = qc.drop(columns=['PC' + str(i + 1) for i, col in enumerate(qc.columns) if i >= pc_max])
    print('qc', qc)
    return qc


def extract_plink_assoc_pcfromukb_(args):
    # actually, eigenvec has been aligned to fam
    #fam = read_fam(args.fplink)
    fam = read_fam(args.fplink, chrom_split=args.chrom_split)
    #eigenvec = read_eigenvec(args.fplink, args.plink_assoc_cov_use_sscore)
    eigenvec = load_qc_ukb(args.fqc_pc)
    print('eigenvec', eigenvec)
    cov = pd.read_csv(args.fcov, sep='\\s+',
                      header=0)
    cov[['fid', 'iid']] = cov[['fid', 'iid']].astype(np.int)

    fam = fam.set_index('iid')
    #eigenvec = eigenvec.set_index('iid')
    cov = cov.set_index('iid')

    #eigenvec = eigenvec.drop(columns=['fid'])

    index_intersection = fam.index.intersection(cov.index)
    cov = cov.loc[index_intersection]

    #assert len(eigenvec) == len(fam), "sth wrong"

    cov_assoc = cov.join(eigenvec, how='inner')

    assert len(cov_assoc) == len(fam), "sth wrong"

    print("fam", fam.head())

    cov_assoc = cov_assoc.loc[fam.index]

    print("cov_assoc", cov_assoc.head())

    cov_assoc = cov_assoc.reset_index()

    col_order = ['fid', 'iid', 'sex', 'age'] + ['PC' + str(i + 1) for i in range(10)]
    cov_assoc = cov_assoc[col_order]

    print(args.fout + '.eigenvec_sex_age.cov')
    write_cov(args.fout + '.eigenvec_sex_age.cov', cov_assoc)


def extract_plink_assoc_cov_(args):
    # actually, eigenvec has been aligned to fam
    fam = read_fam(args.fplink)
    eigenvec = read_eigenvec(args.fplink, args.plink_assoc_cov_use_sscore)
    cov = pd.read_csv(args.fcov, sep='\\s+',
                      header=0)
    cov[['fid', 'iid']] = cov[['fid', 'iid']].astype(np.int)

    fam = fam.set_index('iid')
    eigenvec = eigenvec.set_index('iid')
    cov = cov.set_index('iid')

    eigenvec = eigenvec.drop(columns=['fid'])

    index_intersection = fam.index.intersection(cov.index)
    cov = cov.loc[index_intersection]

    assert len(eigenvec) == len(fam), "sth wrong"

    cov_assoc = cov.join(eigenvec, how='inner')

    assert len(cov_assoc) == len(fam), "sth wrong"

    print("fam", fam.head())

    cov_assoc = cov_assoc.loc[fam.index]

    print("cov_assoc", cov_assoc.head())

    cov_assoc = cov_assoc.reset_index()

    write_cov(args.fplink + '.eigenvec_sex_age.cov', cov_assoc)


def standardize_cov(args):
    cov = pd.read_csv(args.fcov, sep='\\s+',
                      header=0)
    cov.loc[:, ['fid', 'iid']] = cov[['fid', 'iid']].astype(np.int)
    #cov = cov.set_index('iid')
    cols = list(cov.columns)
    cols.remove('fid')
    cols.remove('iid')

    cov.loc[:, cols] = (cov[cols] - cov[cols].mean()) / cov[cols].std()

    write_cov(args.fcov + '.standardize.cov', cov)


def change_order_cov_(args):
    # samples of cov are aligned to ukb_imp_one.fam
    # change the order to ukb_imp_chr%.fam

    # for cv0
    fam = read_fam(args.fplink,
                   chrom_split=args.chrom_split)

    # fam = read_fam(args.fplink+'.obs.cv0.tr.unrelated',
    #               chrom_split=args.chrom_split)

    cov = pd.read_csv(args.fcov, sep='\\s+',
                      header=0)
    cov[['fid', 'iid']] = cov[['fid', 'iid']].astype(np.int)

    fam = fam.set_index('iid')
    cov = cov.set_index('iid')

    print("cov", cov.head())
    print("fam", fam.head())

    cov = cov.loc[fam.index, :]

    print("aligned cov", cov.head())

    cov = cov.reset_index()

    print("write to: ", args.fout)
    # indicate fout_com+'...' as fout
    write_cov(args.fout, cov)
    #write_cov(args.fout+'.obs.cv0.tr.unrelated.eigenvec_sex_age.cov', cov)


def change_order_cov_com2one(args):
    # samples of cov are aligned to ukb_imp_one.fam
    # change the order to ukb_imp_chr%.fam

    # for plink one

    for i_cv in range(args.n_cv):
        fam = read_fam(args.fplink + '.obs.cv' + str(i_cv) + '.tr.unrelated',
                       chrom_split=args.chrom_split)

        cov = pd.read_csv(args.fcov + '.obs.cv' + str(i_cv) + '.tr.unrelated.eigenvec_sex_age.cov', sep='\\s+',
                          header=0)
        cov[['fid', 'iid']] = cov[['fid', 'iid']].astype(np.int)
        cov_cols = cov.columns

        fam = fam.set_index('iid')
        cov = cov.set_index('iid')

        print("cov", cov.head())
        print("fam", fam.head())

        cov = cov.loc[fam.index, :]

        print("aligned cov", cov.head())

        cov = cov.reset_index()

        cov = cov.loc[:, cov_cols]

        print("write to: ", args.fout)
        # indicate fout_com+'...' as fout
        write_cov(args.fout + '.obs.cv' + str(i_cv) + '.tr.unrelated.eigenvec_sex_age.cov', cov)
        #write_cov(args.fout+'.obs.cv0.tr.unrelated.eigenvec_sex_age.cov', cov)

    # add test if needed


def add_a2_plinkassoc_(args):

    for chrom in range(1, 23):
        # load bim
        print("chrom", chrom)
        bim = read_bim(args.fplink, chrom_split=True, chrom=chrom)

        print("bim.head()", bim.head())

        fassoc = args.fplinkassoc.replace('%', str(chrom))
        assoc = pd.read_csv(fassoc, sep='\\s+', header=0)
        print("assoc.head()", assoc.head())
        assoc['A2'] = 'N'

        assert len(bim) == len(assoc), "sth wrong"

        print("original bim", len(bim))

        filt = (assoc['A1'] == bim['A1'])
        assoc.loc[filt, 'A2'] = bim.loc[filt, 'A2']

        print("assoc.A1==bim.A1", len(assoc[filt]))

        filt = (assoc['A1'] == bim['A2'])
        assoc.loc[filt, 'A2'] = bim.loc[filt, 'A1']

        print("assoc.A1==bim.A2", len(assoc[filt]))

        assert len(assoc['A2'] == 'N') > 0, "sth wrong"

        print("bfr exclude p is nan", len(assoc))

        assoc = assoc[~assoc['P'].isnull()]
        print("afr exclude p is nan", len(assoc))

        print("save to ", fassoc + '.A2')
        assoc.to_csv(fassoc + '.A2', sep='\t', index=False)


def adjust_pc_projection_sscore(args):
    """
    use eigeval, eigenvec, sscore_training to derive relationship between eigenvec and sscore.
    sscore_training is --score on the same training dataset of eigenvec.
    eigenvec is used as covariant in adaboost.
    """
    eigenval = pd.read_csv(args.feigenval, header=None)
    print("eigenval", eigenval)
    eigenvec = pd.read_csv(args.feigenvec, sep='\t')
    print("eigenvec", eigenvec)
    sscore_training = pd.read_csv(args.fsscore_training, sep='\t')
    print("sscore_training", sscore_training)
    sscore = pd.read_csv(args.fsscore, sep='\t')
    # print("sscore",sscore)

    print("columns", eigenvec.columns)
    sscore.columns = eigenvec.columns

    eigenval_np = eigenval.values
    eigenvec_np = eigenvec.iloc[:, 2:].values
    sscore_training_np = sscore_training.iloc[:, 2:].values

    # ver.1
    print("ver1,vec=sscore/np.sqrt(eigenval)")
    eigenvec_val_ = sscore_training_np / np.sqrt(eigenval_np.T)
    # eigenvec_val_=np.sqrt(eigenval_np.T)*eigenvec_np
    print("eigenvec_val_", eigenvec_val_)

    print("np", eigenvec_val_)
    print("np", eigenvec_np)
    if np.allclose(eigenvec_val_, eigenvec_np, rtol=1e-1):
        print("ver.1 is right")
        sscore.iloc[:, 2:] = sscore.iloc[:, 2:] / np.sqrt(eigenval_np.T)
        print("sscore", sscore)
        sscore.to_csv(args.fout, sep='\t', index=False)

        return

    # ver.2
    print("ver2,vec=-sscore/np.sqrt(eigenval)")
    eigenvec_val_ = -sscore_training_np / np.sqrt(eigenval_np.T)
    # eigenvec_val_=-np.sqrt(eigenval_np.T)*eigenvec_np
    print("eigenvec_val_", eigenvec_val_)

    # if np.allclose(eigenvec_val_, sscore_training_np,atol=1e-3):
    if np.allclose(eigenvec_val_, eigenvec_np, rtol=1e-1):
        print("ver.2 is right")
        sscore.iloc[:, 2:] = -sscore.iloc[:, 2:] / np.sqrt(eigenval_np.T)
        print("sscore", sscore)
        sscore.to_csv(args.fout, sep='\t', index=False)

        return

    print("no ver. fit")

    # sscore_eigen


def main():
    args = argument()

    # extract samples sent to out e-mail
    # prepare withdrawl_iid.txt
    # and remove

    if args.valid:
        extract_samples_(args)

    if args.cov:
        write_cov_(args)

    if args.case:
        write_case_(args)

    if args.case_multi:
        write_case_multi_(args)

    if args.snps:
        extract_snps_(args)

    if args.snps_dup_impinfo:
        extract_snps_dup_impinfo_(args)

    if args.snps_hm3:
        extract_snps_hm3(args)

    if args.assoc_hm3:
        extract_assoc_hm3(args)

    if args.snps_hm3_freq:
        extract_snps_hm3_freq(args)

    if args.snps_hm3_sida_rs:
        extract_snps_hm3_sida_rs(args)

    if args.sida_rs:
        map_sida_rs(args)

    if args.gene_exon:
        create_gene_exon(args)

    if args.samples_female:
        extract_samples_female_(args)

    if args.unrelated:
        extract_samples_unrelated_(args)

    if args.unrelated_va:
        extract_samples_unrelated_va_(args)

    if args.unrelated_va12:
        extract_samples_unrelated_va12_(args)

    if args.unrelated_test:
        extract_samples_unrelated_test_(args)

    if args.training_propcase:
        extract_training_propcase_(args)

    if args.plink_assoc_cov_pcfromukb:
        extract_plink_assoc_pcfromukb_(args)

    if args.split_covs:
        split_covs(args)

    if args.align_covs:
        align_covs(args)

    if args.align_multiphe:
        align_multiphe(args)

    if args.zeroone_multiphe:
        # for snpnet
        zeroone_multiphe(args)

    if args.plink_assoc_cov_pcfromukb_com2one:
        change_order_cov_com2one(args)

    if args.plink_assoc_cov:
        extract_plink_assoc_cov_(args)

    if args.plink_assoc_cov_one2chrom:
        change_order_cov_(args)

    if args.add_a2_plinkassoc:
        add_a2_plinkassoc_(args)

    if args.adjust_pc_projection:
        adjust_pc_projection_sscore(args)

    if args.down_impinfo:
        extract_snps_down_impinfo(args)

    if args.down_cal:
        extract_snps_down_cal(args)

    if args.cov_standardize:
        standardize_cov(args)

    print('Done!')


if __name__ == '__main__':
    main()
