
"""
simulation
"""


import argparse
import pathlib
import os
import hashlib
import pandas as pd
import numpy as np
from logging import getLogger
from ...system.logger import logger_setting

from ..io_genot import plink
from ..snv import snvio as snvio
# from ..snv import snv as snvop


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--create-wgt',
                        action='store_true',
                        help='')

    parser.add_argument('--create-phe',
                        action='store_true',
                        help='')

    parser.add_argument('--dout',
                        default=None,
                        help='output file')

    parser.add_argument('--fgenot',
                        default=None,
                        help='input prefix for genotype file')

    parser.add_argument('--genot-format',
                        default='plink2vzs',
                        help='plink2, plink, etc.')

    parser.add_argument('--ffreq',
                        default=None,
                        help='')


    parser.add_argument('--liability-scale',
                        action='store_true',
                        help='')


    parser.add_argument('--randomn',
                        help='')

    parser.add_argument('--seed',
                        type=int,
                        default=None,
                        help='seed of random_state')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def create_wgt_simu_random(snv, h2_add, h2_dom, ncausal,liability_scale,prevalence):
    '''
    causal variants position are random
    '''
    nsnv = len(snv)

    perm = np.random.permutation(nsnv)
    use_snv_index = perm[:ncausal]

    filt_non0 = np.zeros(nsnv, dtype=bool)
    filt_non0[use_snv_index] = True

    # filt_non0 = np.random.rand(nsnv)
    # prop = ncausal/nsnv
    # filt_non0 = (filt_non0 < prop)
    # nnon0 = np.count_nonzero(filt_non0)
    # print('nnon0/nsnv:',nnon0,nsnv)
    # if nnon0==0:
    # raise RuntimeError

    # freq_non0 = snv['alt_frq'].to_numpy()[filt_non0]
    # #snv * 3
    # dosage_score=freq_non0
    wgt = snv.loc[filt_non0, ['id', 'chrom', 'pos', 'ref', 'alt', 'alt_frq']].copy().reset_index(drop=True)

    # add_std = np.sqrt(h2_add / (nsnv * prop))
    add_std = np.sqrt(h2_add / ncausal)
    add_effects_normalized = np.random.normal(loc=0.0, scale=add_std, size=ncausal)

    # dom_std = np.sqrt(h2_dom / (nsnv * prop))
    dom_std = np.sqrt(h2_dom / ncausal)
    dom_effects_normalized = np.random.normal(loc=0.0, scale=dom_std, size=ncausal)

    wgt.loc[:, 'add'] = add_effects_normalized
    # wgt.loc[filt_non0, 'add'] = add_effects_normalized
    wgt.loc[:, 'dom'] = dom_effects_normalized
    # wgt.loc[:, 'alt_frq'] = snv.loc[filt_non0, 'alt_frq']

    alt_frq = wgt['alt_frq']
    # print('alt_frq', alt_frq)
    ref_frq = 1.0 - wgt['alt_frq']
    add_coef = wgt['add'] * (1 / np.sqrt(2 * alt_frq * (1 - alt_frq)))
    wgt.loc[:, 'score0'] = wgt['dom'] * (-alt_frq / ref_frq)
    wgt.loc[:, 'score1'] = add_coef * 1 + wgt['dom']
    wgt.loc[:, 'score2'] = add_coef * 2 + wgt['dom'] * (-ref_frq / alt_frq)

    wgt = wgt.rename(columns={'alt': 'a1', 'ref': 'a2'})

    return wgt


def hash_from_string(s):
    # 8 digits should be fine
    return int(hashlib.sha256(s.encode('utf-8')).hexdigest(), 16) % 10**8


def run_simu_wgt(dout, fgenot, genot_format, ffreq,liability_scale, randomn, seed):
    # no reproducability
    # np.random.seed(seed)

    snv = plink.load_snv(fgenot, genot_format)

    freq = plink.load_freq_plink2(ffreq)
    logger.debug('freq: {}'.format(freq.info()))
    freq_check = freq[['id', 'ref', 'alt']].set_index('id')
    snv_check = snv.set_index(('id'))
    if (freq_check['ref'] != snv_check['ref']).any() or (freq_check['alt'] != snv_check['alt']).any():
        raise RuntimeError('unmatch allele')

    freq = freq[['id', 'alt_frq']]
    snv = pd.merge(snv, freq, on='id', how='left')

    #h2s = [('0.5', '0.0'),
    #       ('0.4', '0.1'),
    #       ('0.45', '0.05'),
    #       ('0.0', '0.5'),
    #       ('0.1', '0.0'),
    #       ('0.08', '0.02'),
    #       ('0.09', '0.01'),
    #       ('0.0', '0.1'),
    #       ('0.2', '0.0'),
    #       ('0.16', '0.04'),
    #       ('0.18', '0.02'),
    #       ('0.0', '0.2'),
    #       ('0.25', '0.0'),
    #       ('0.2', '0.05'),
    #       ('0.0', '0.25')
    #       ]

    h2s = [
           ('0.1', '0.0'),
           ('0.08', '0.02'),
           ('0.09', '0.01'),
           ('0.0', '0.1'),
           ('0.25', '0.0'),
           ('0.2', '0.05'),
           ('0.0', '0.25')
           ]


    #  100/1M=1e-4
    # 1000/1M=1e-3
    #ncausals = ['100', '1000']
    ncausals = ['10000']
    for ncausal in ncausals:
        for (h2_add, h2_dom) in h2s:
            logger.debug('para: h2add, h2dom, ncausal: {}, {}, {}'.format(h2_add, h2_dom, ncausal))
            for random_i in range(int(randomn)):

                model = 'h2add-' + h2_add + '_h2dom-' + h2_dom + '_ncausal-' + ncausal

                # para_tuple = (h2_add, h2_dom, ncausal)
                seed_para = (seed + hash_from_string(model) + random_i) % 2**32
                # NEVER USE hash() since the value changes when rerun python
                # seed_para = (seed + hash(para_tuple) + random_i) % (2**32)
                logger.debug('seed_para: {}'.format(seed_para))
                np.random.seed(seed_para)

                # random_0 for plink column name; random-0 raise error
                fwgt = dout / model / 'true-wgt' / ('random_' + str(random_i) + '.wgt')
                #fwgt = dout / model / 'true-wgt' / ('random-' + str(random_i) + '.wgt')
                # fwgt = dout / ('h2add-' + h2_add + '_h2dom-' + h2_dom + '_ncausal-' + ncausal) / 'true-wgt' / ('random-' + str(random_i) + '.wgt')
                wgt = create_wgt_simu_random(snv, float(h2_add), float(h2_dom), int(ncausal),liability_scale,)
                logger.debug('wgt: {}'.format(wgt))
                os.makedirs(fwgt.parent, exist_ok=True)
                wgt.to_csv(fwgt, sep='\t', index=False)


def parse_model(model):
    ps = model.split('_')
    if len(ps) != 3:
        raise RuntimeError

    h2add = ps[0].split('-')[1]
    h2dom = ps[1].split('-')[1]
    ncausal = ps[2].split('-')[1]

    return (h2add, h2dom, ncausal)


def create_score(score_genot, h2_add, h2_dom):
    logger.debug('h2add, h2dom {}, {}'.format(h2_add, h2_dom))

    var = np.var(score_genot.loc[:, 'score'])
    # check if var ~ h2add+h2dom
    logger.debug('h2, genot var: {} vs {:.4}'.format(h2_add + h2_dom, var))

    # use on create phenotype
    h2_remain = 1.0 - (h2_add + h2_dom)

    sample_n = len(score_genot)

    sample_score = np.random.normal(loc=0.0, scale=np.sqrt(h2_remain), size=sample_n)

    score_phe = score_genot.copy()
    score_phe.loc[:, 'liability'] = score_genot.loc[:, 'score'] + sample_score

    logger.debug('env var: {:.4}'.format(np.var(sample_score)))
    logger.debug('liability var: {:.4}'.format(np.var(score_phe['liability'])))

    prevalence = 0.1

    score_phe_sort = np.sort(score_phe['liability'])[::-1]
    score_thre = score_phe_sort[int(sample_n * prevalence)]

    score_phe.loc[:, 'status'] = (score_phe['liability'] > score_thre).astype(int)

    # print('ncase / sample_n', np.count_nonzero(score_phe.loc[:, 'status'] == 1), sample_n)

    return score_phe


def run_simu_phe(dout, seed):
    # no replicatability
    # np.random.seed(seed)

    ds = os.listdir(dout)

    for d in ds:
        model = d
        para_tuple = parse_model(model)
        # (h2_add, h2_dom, _) = parse_model(model)
        (h2_add, h2_dom, _) = para_tuple

        logger.debug('\n')
        logger.debug('model: {}'.format(model))

        dscore = dout / d / 'score-genot'
        for frandom_name in os.listdir(dscore):
            if not frandom_name.endswith('.score'):
                continue

            # logger.debug('framdom_name: {}'.format(frandom_name))
            score_genot = pd.read_csv(dscore / frandom_name, sep='\t', header=0)

            random_name = pathlib.Path(frandom_name).stem

            random_i = int(random_name.split('_')[1])
            seed_para = (seed + hash_from_string(model) + random_i) % 2**32
            logger.debug('seed_para: {}'.format(seed_para))
            np.random.seed(seed_para)

            pheno = create_score(score_genot, float(h2_add), float(h2_dom))
            # print('pheno', pheno)

            fpheno = dout / d / 'phe' / 'raw' / (random_name + '.phe')
            os.makedirs(fpheno.parent, exist_ok=True)
            pheno.to_csv(fpheno, sep='\t', index=False)

    # concat
    for d in ds:
        draw = dout / d / 'phe' / 'raw'
        frandom_names = os.listdir(draw)
        phenos = [pd.read_csv(draw / fphe, sep='\t') for fphe in frandom_names]

        pheno_concat = phenos[0].loc[:, ['id']].copy()

        for (pheno, frandom_name) in zip(phenos, frandom_names):
            random_name = frandom_name.split('.')[0]
            # plink2 does not allow '-'
            random_name = random_name.replace('-', '_')
            pheno_concat = pd.merge(pheno_concat, pheno[['id', 'status']].rename(columns={'status': random_name}), on='id')
            print('pheno_concat', pheno_concat)

        pheno_concat = pheno_concat.rename(columns={'id': 'IID'})
        fpheno_concat = dout / d / 'phe' / 'phe.phe'
        pheno_concat.to_csv(fpheno_concat, sep='\t', index=False)


def main():
    args = argument()

    if args.create_wgt:
        run_simu_wgt(pathlib.Path(args.dout),
                     args.fgenot, args.genot_format,
                     args.ffreq,
                     args.liability_scale,
                     args.randomn, args.seed)

    if args.create_phe:
        run_simu_phe(pathlib.Path(args.dout), args.seed)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
