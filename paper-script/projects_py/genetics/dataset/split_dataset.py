"""
For
- QC on samples or snps.
- split dataset.

input is plink file .fam or .bim.
output is extracted samples file or snps files: "${fout}.samples", "${fout}.snps"

After running this program, you would
`$ plink --bfile $fplink --extract ${fplink}.samples --keep ${fplink}.snps --out ${fplink_out}`

TODO: use df to output $fsamples since last \n will be added
"""

import os
import time
import argparse

import random

import numpy as np
import pandas as pd

from sklearn.model_selection import ShuffleSplit,StratifiedShuffleSplit, KFold,StratifiedKFold

# import io_gwas as iog
# import validate_ukb


def argument():
    parser = argparse.ArgumentParser(
        prog='gwas',
        usage='gwas',
        description='description',
        epilog='description end',
        add_help=True,
    )


    parser.add_argument('--multiphe',
                        action='store_true',
                        help="split datasets common for phes")

    parser.add_argument('--samples',
                        action='store_true',
                        help="qc on samples")

    parser.add_argument('--cov',
                        action='store_true',
                        help="split cov")

    parser.add_argument('--cov_eigenvec_sex_age',
                        action='store_true',
                        help="split cov")

    parser.add_argument('--unrelated',
                        action='store_true',
                        help="**TMP** only for some of them; like --split_training_12")

    parser.add_argument(
        '--fplink',
        help="prefix of plink file without '.bed', if $chrom_split, indicate '%' as chromosome number e.g. 'ukb_imp_chr%_test'.")

    parser.add_argument('--fcov',
                        default=None,
                        help="cov file")

    parser.add_argument('--chrom_split',
                        action='store_true',
                        help="chromosomes are split in fplink")
    """
    # use fam
    parser.add_argument('--fcase',
                        help="file name of list of case samples with one per line. \n\
 also accept space / tab - delimited text file with family IDs in the first column and within - family IDs in the second column, which is the same format as input of plink - -keep command ")
    """

    parser.add_argument(
        '--split_observed_test',
        action='store_true',
        help="split samples into observed dataset (=training dataset and validation dataset) and test dataset.")

    parser.add_argument(
        '--split_training_validation',
        action='store_true',
        help="split observed datset into training dataset and validation dataset if with --split_observed_test. if not, split whole dataset into them.")

    parser.add_argument('--split_observed',
                        action='store_true',
                        help="split samples into observed dataset. only for .cov")

    parser.add_argument('--split_validation_12',
                        action='store_true',
                        help="split validation samples into 1,2")

    parser.add_argument('--split_training_12',
                        action='store_true',
                        help="split trainig samples into 1,2")

    parser.add_argument('--prop_test',
                        type=float,
                        default=0.2,
                        help="propotion of test to whole dataset")

    parser.add_argument('--prop_12',
                        type=float,
                        default=0.1,
                        help="use for validation or test")

    parser.add_argument('--n_cv',
                        type=int,
                        default=5,
                        help='numbers to split for cross-validation')

    parser.add_argument('--seed',
                        type=int,
                        default=None,
                        help='seed of random_state')
    """
    parser.add_argument('--fout',
                        default=None,
                        help="output prefix for each chromosome.")
                        """

    parser.add_argument('--fout_com',
                        default=None,
                        help="output prefix common for chromosome.")

    """
    parser.add_argument('--fout',
                        help="output prefix. \
a\
output ${out}.samples and ${out}.snps. if split_training_validation, output ${out}.train.samples and ${out}.vali.samples and ${out}.snps")
    """

    args = parser.parse_args()

    print(args)

    return args


def read_fam(fplink, chrom_split=False):
    cols = 'fid,iid,N1,N2,sex,phe'.split(',')
    if chrom_split:
        fplink = fplink.replace('%', str(1))
    pheno = pd.read_csv(fplink + '.fam', sep='\\s+',
                      header=None,
                      names=cols)
    # check if not 0 or 1
    pheno['phe'] = pheno['phe'] - 1
    return pheno


def read_case(fcase):
    case = pd.read_csv(fcase, sep='\\s+',
                       header=None)
    return case


def add_case(pheno, case):
    pheno['phe'] = 0

    if len(case.columns) == 2:
        cols = ['fid', 'iid']
    elif len(case.columns) == 1:
        cols = ['iid']
    else:
        raise IOError("number of columns is not 1 or 2")

    case.columns = cols
    case = case.set_index(cols)
    pheno = pheno.set_index(cols)
    pheno.loc[case.index, 'phe'] = 1
    pheno = pheno.reset_index(drop=False)

    return pheno


def split_samples_multiphe(phe, prop=0.2, seed=51):
    pheno_vals = phe['phe'].values

    sss = ShuffleSplit(
        n_splits=1, test_size=prop, random_state=seed)
    indexs = list(sss.split(np.zeros_like(pheno_vals)))[0]
    indexs = [np.sort(v) for v in indexs]
    print("indexs[0]", indexs[0][:10])
    print("indexs[1]", indexs[1][:10])
    return indexs

def split_samples_multiphe_cv(pheno, n_cv=5, seed=52):
    pheno_vals = pheno['phe'].values
    print("pheno", pheno_vals)

    skf = KFold(n_splits=n_cv, shuffle=True, random_state=seed)
    # list of n_cv tuples, each contains (index_tr),(index_va)
    indexss = list(skf.split(np.zeros_like(pheno_vals), pheno_vals))
    indexss = [[np.sort(v) for v in vs] for vs in indexss]
    print("indexss[0][0]", indexss[0][0][:10])
    return indexss



def split_samples(pheno, prop=0.2, seed=51):
    pheno_vals = pheno['phe'].values

    sss = StratifiedShuffleSplit(
        n_splits=1, test_size=prop, random_state=seed)
    indexs = list(sss.split(np.zeros_like(pheno_vals), pheno_vals))[0]
    indexs = [np.sort(v) for v in indexs]
    print("indexs[0]", indexs[0][:10])
    print("indexs[1]", indexs[1][:10])
    return indexs


def split_samples_cv(pheno, n_cv=5, seed=52):
    pheno_vals = pheno['phe'].values
    print("pheno", pheno_vals)

    skf = StratifiedKFold(n_splits=n_cv, shuffle=True, random_state=seed)
    # list of n_cv tuples, each contains (index_tr),(index_va)
    indexss = list(skf.split(np.zeros_like(pheno_vals), pheno_vals))
    indexss = [[np.sort(v) for v in vs] for vs in indexss]
    print("indexss[0][0]", indexss[0][0][:10])
    return indexss


def write_index_direct(fout, pheno):

    print("phe", len(pheno))
    #iid_ = list(phe.iloc[index]['iid'].values)
    iid_ = list(pheno['iid'].values)
    print("iid_", len(iid_))
    print("iid_", iid_[:5])
    iid_ = [str(v) + '\t' + str(v) for v in iid_]
    iid_ = '\n'.join(iid_)
    with open(fout, 'w') as f:
        f.write(iid_)

        
def write_index_pd(fout, pheno, indexs, names=['obs', 'test'], i_cv=None):
    # this way adds '\n' to the last

    def write_index_(fout, index):
        print("phe", len(pheno))
        print("index", len(index))

        pheno.iloc[index][['fid','iid']].to_csv(fout,header=None,index=False,sep='\t')
        
        #iid_ = list(phe.iloc[index]['iid'].values)
        #print("iid_", len(iid_))
        #print("iid_", iid_[:5])
        #iid_ = [str(v) + '\t' + str(v) for v in iid_]
        #iid_ = '\n'.join(iid_)
        #with open(fout, 'w') as f:
        #    f.write(iid_)

    if i_cv is None:
        """
        fout_1 = fout + '.' + names[0] + '.samples'
        fout_2 = fout + '.' + names[1] + '.samples'
        """
        fouts_ = [fout + '.' + name + '.samples' for name in names]
    else:
        """
        fout_1 = fout + '.cv' + str(i_cv) + '.' + names[0] + '.samples'
        fout_2 = fout + '.cv' + str(i_cv) + '.' + names[1] + '.samples'
        """
        fouts_ = [fout + '.cv' + str(i_cv) + '.' +
                  name + '.samples' for name in names]

    for fout_, index in zip(fouts_, indexs):
        print("output to:", fout_)
        print("len(indexs)", len(indexs))
        print("len(index)", len(index))
        write_index_(fout_, index)




def write_index(fout, pheno, indexs, names=['obs', 'test'], i_cv=None):

    def write_index_(fout, index):
        print("phe", len(pheno))
        print("index", len(index))
        iid_ = list(pheno.iloc[index]['iid'].values)
        print("iid_", len(iid_))
        print("iid_", iid_[:5])
        iid_ = [str(v) + '\t' + str(v) for v in iid_]
        iid_ = '\n'.join(iid_)
        with open(fout, 'w') as f:
            f.write(iid_)

    if i_cv is None:
        """
        fout_1 = fout + '.' + names[0] + '.samples'
        fout_2 = fout + '.' + names[1] + '.samples'
        """
        fouts_ = [fout + '.' + name + '.samples' for name in names]
    else:
        """
        fout_1 = fout + '.cv' + str(i_cv) + '.' + names[0] + '.samples'
        fout_2 = fout + '.cv' + str(i_cv) + '.' + names[1] + '.samples'
        """
        fouts_ = [fout + '.cv' + str(i_cv) + '.' +
                  name + '.samples' for name in names]

    for fout_, index in zip(fouts_, indexs):
        print("output to:", fout_)
        print("len(indexs)", len(indexs))
        print("len(index)", len(index))
        write_index_(fout_, index)

    """
    # use com
    if not chrom_split:
        for fout_, index in zip(fouts_, indexs):
            write_index_(fout_, index)
    else:
        for fout_, index in zip(fouts_, indexs):
            for chrom in range(1, 23):
                fout_ = fout_.replace('%', str(chrom))
                write_index_(fout_, index)
                """


def read_index(fin, names=['obs', 'test'], i_cv=None):

    def read_index_(fin):
        df_index = pd.read_csv(
            fin, sep='\t', header=None, names=['fid', 'iid'])
        index = df_index['iid'].values
        return index

    if i_cv is None:
        fouts_ = [fout + '.' + name + '.samples' for name in names]
    else:
        fouts_ = [fout + '.cv' + str(i_cv) + '.' +
                  name + '.samples' for name in names]

    indexs = [read_index_(fin) for fin in fouts_]
    return indexs


def read_cov(fin):
    df = pd.read_csv(fin, sep='\t', header=0)
    return df


def write_cov(fout, cov):
    cov.to_csv(fout, sep='\t',
               index=False)


"""
def write_split_cov(fout, cov, indexs, names=['obs', 'test']):
    def write_cov_(fout, cov_):
        pass

    if i_cv is None:
        fouts_ = [fout + '.' + name + '.samples' for name in names]
    else:
        fouts_ = [fout + '.cv' + str(i_cv) + '.' +
                  name + '.samples' for name in names]

    for fout_, index in zip(fouts_, indexs):
        # extract cov_ with index
        print("output to:", fout_)
        write_cov_(fout_, cov_)
        """


def extract_cov_write(fout, cov, fplink, chrom_split=False):
    print("fplink: ", fplink)
    fam = read_fam(fplink, chrom_split)
    print("fam", fam.head())
    cov_ = cov.loc[fam['iid']]
    cov_ = cov_.reset_index()
    print("cov_ should align to fam", cov_.head())
    write_cov(fout, cov_)


def main():
    args = argument()

    if args.samples:
        phe = read_fam(args.fplink, args.chrom_split)
        print("phe", phe.head())
        # phe = pd.read_csv(args.fplink)
        # phenotype is already in phe
        # case = read_case(args.fcase)
        # phe = add_case(phe, case)

        if args.multiphe:
            if args.split_observed_test:
                # first, split whole dataset into observed and test
                indexs_obs_ts = split_samples_multiphe(
                    phe, prop=args.prop_test, seed=args.seed)
                write_index_pd(args.fout_com, phe, indexs_obs_ts,
                            names=['obs', 'test'])

                if args.split_training_validation:
                    # next, split only observed datset into training and validation
                    phe_obs = phe.iloc[indexs_obs_ts[0]]
                    print("phe_obs", phe_obs.head())
                    phe_obs = phe_obs.reset_index(drop=True)
                    print("resetindex", phe_obs.head())
                    indexss_cv = split_samples_multiphe_cv(
                        phe_obs, n_cv=args.n_cv, seed=args.seed + 1)
                    for i_cv, indexs in enumerate(indexss_cv):
                        write_index_pd(args.fout_com + '.obs', phe_obs, indexs,
                                    names=['tr', 'va'], i_cv=i_cv)


        else:

            if args.split_observed_test:
                # first, split whole dataset into observed and test
                indexs_obs_ts = split_samples(
                    phe, prop=args.prop_test, seed=args.seed)
                write_index(args.fout_com, phe, indexs_obs_ts,
                            names=['obs', 'test'])

                if args.split_training_validation:
                    # next, split only observed datset into training and validation
                    phe_obs = phe.iloc[indexs_obs_ts[0]]
                    print("phe_obs", phe_obs.head())
                    phe_obs = phe_obs.reset_index(drop=True)
                    print("resetindex", phe_obs.head())
                    indexss_cv = split_samples_cv(
                        phe_obs, n_cv=args.n_cv, seed=args.seed + 1)
                    for i_cv, indexs in enumerate(indexss_cv):
                        write_index(args.fout_com + '.obs', phe_obs, indexs,
                                    names=['tr', 'va'], i_cv=i_cv)

            # do not use
            #elif args.split_training_validation:
            #    # split whole dataset into training and validation
            #    indexss_cv = split_samples_cv(
            #        phe, n_cv=args.n_cv, seed=args.seed + 1)
            #    for i_cv, indexs in enumerate(indexss_cv):
            #        write_index(args.fout_com, phe, indexs,
            #                    names=['tr', 'va'], i_cv=i_cv)

            # elif args.split_training_12:
                # args.fout_com+'.obs.cv'+str(i_cv)+'.tr.unrelated.samples'
            elif args.split_validation_12:

                def read_samples(fsample):
                    cols = 'fid,iid'.split(',')
                    sample = pd.read_csv(fsample, sep='\\s+',
                                        header=None,
                                        names=cols)
                    return sample

                for i_cv in range(args.n_cv):
                    print("i_cv", i_cv)
                    phe = read_fam(args.fplink, args.chrom_split)
                    print("phe", phe.head())
                    fsample = args.fout_com + '.obs.cv' + str(i_cv) + '.va.samples'
                    samples_va = read_samples(fsample).loc[:, 'iid']
                    phe_va = phe.set_index('iid').loc[samples_va].reset_index(drop=False)
                    # phe_va=phe.loc[phe['iid'].isin(samples_va),:].reset_index(drop=True)
                    print("phe_va", phe_va.head())

                    indexs_va12 = split_samples(
                        phe_va, prop=args.prop_12, seed=args.seed + 2)

                    # for va1, va2
                    for va_name, indexs_va in zip(['1', '2'], indexs_va12):
                        phe_va_write = phe_va.iloc[indexs_va]

                        fout_write = args.fout_com + '.obs.cv' + str(i_cv) + '.va' + va_name + '.samples'
                        print("fout_write", fout_write)
                        write_index_direct(fout_write, phe_va_write)

            elif args.split_training_12:

                def read_samples(fsample):
                    cols = 'fid,iid'.split(',')
                    sample = pd.read_csv(fsample, sep='\\s+',
                                        header=None,
                                        names=cols)
                    return sample

                for i_cv in range(args.n_cv):
                    print("i_cv", i_cv)
                    phe = read_fam(args.fplink, args.chrom_split)
                    print("phe", phe.head())
                    if args.unrelated:
                        fsample = args.fout_com + '.obs.cv' + str(i_cv) + '.tr.unrelated.samples'
                    else:
                        fsample = args.fout_com + '.obs.cv' + str(i_cv) + '.tr.samples'
                    samples_va = read_samples(fsample).loc[:, 'iid']
                    phe_va = phe.set_index('iid').loc[samples_va].reset_index(drop=False)
                    # phe_va=phe.loc[phe['iid'].isin(samples_va),:].reset_index(drop=True)
                    print("phe_va", phe_va.head())

                    indexs_va12 = split_samples(
                        phe_va, prop=args.prop_12, seed=args.seed + 3)

                    # for va1, va2
                    for va_name, indexs_va in zip(['1', '2'], indexs_va12):
                        phe_va_write = phe_va.iloc[indexs_va]

                        if args.unrelated:
                            fout_write = args.fout_com + '.obs.cv' + str(i_cv) + '.tr' + va_name + '.unrelated.samples'
                        else:
                            fout_write = args.fout_com + '.obs.cv' + str(i_cv) + '.tr' + va_name + '.samples'
                        print("fout_write", fout_write)
                        write_index_direct(fout_write, phe_va_write)

    if args.cov:
        # order is sorted to .fam
        cov = read_cov(args.fcov)
        cov = cov.set_index('iid')

        # is this ok??
        if args.split_observed_test:
            # indexs = read_index(args.fout_com)
            for name in ['obs', 'test']:
                print("name", name)
                extract_cov_write(args.fout_com + '.' + name + '.cov',
                                  cov=cov,
                                  fplink=args.fplink + '.' + name,
                                  chrom_split=args.chrom_split)

        if args.split_observed:
            # extract cov in fplink.obs
            # assume for independent samples
            for name in ['obs']:
                extract_cov_write(args.fout_com + '.' + name + '.unrelated' + '.cov',
                                  cov=cov,
                                  fplink=args.fplink + '.' + name,  # TODO: dont you need +'.unrelated'
                                  chrom_split=args.chrom_split)

    if args.cov_eigenvec_sex_age:
        # order is sorted to .fam
        cov = read_cov(args.fcov)
        cov = cov.set_index('iid')

        if args.split_training_12:
            # indexs = read_index(args.fout_com)
            if args.unrelated:
                for i_cv in range(args.n_cv):
                    print("i_cv", i_cv)
                    # for name in ['tr1', 'tr2']:
                    for name in ['tr1']:
                        print("name", name)
                        extract_cov_write(
                            args.fout_com + '.obs.cv' + str(i_cv) + '.' + name +
                            '.unrelated.eigenvec_sex_age.cov',
                            cov=cov,
                            fplink=args.fplink + '.obs.cv' + str(i_cv) + '.' + name + '.unrelated',
                            chrom_split=args.chrom_split)

            """
            # do not create cov file for each cv
            if args.split_training_validation:
                for i_cv in range(n_cv):
                    for name in ['tr', 'va']:
                        extract_cov_write(fplink+'.obs'+str(i_cv) +
                                          '.' + name+'.cov',
                                          cov=cov_,
                                          fplink=fplink+'.obs'+str(i_cv) +
                                          '.' + name,
                                          chrom_split=args.chrom_split)
                                          """

            """
            # first, split whole dataset into observed and test
            indexs_obs_ts = split_samples(
                phe, prop=args.prop_test, seed=args.seed)
            write_index(args.fout_com, indexs_obs_ts,
                        names=['obs', 'test'])

            if args.split_training_validation:
                # next, split only observed datset into training and validation
                phe_obs = phe.iloc[indexs_obs_ts[0]]
                phe_obs = phe_obs.reset_index(drop=True)
                indexss_cv = split_samples_cv(
                    phe_obs, n_cv=args.n_cv, seed=args.seed + 1)
                for i_cv, indexs in enumerate(indexss_cv):
                    write_index(args.fout_com + '.obs', indexs,
                                names=['tr', 'va'], i_cv=i_cv)
                                """

    print('Done!')


if __name__ == '__main__':
    main()
