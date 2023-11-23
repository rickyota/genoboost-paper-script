
import argparse
import os
import pathlib
import yaml
import functools
import numpy as np
import pandas as pd
from logging import getLogger

from ...system.logger import logger_setting
from ...genetics.io_genot import sample as sampleio
from ...genetics.io_genot import score as scoreio
from ...genetics.sample import sample as sampleop
from ..io import filename as iof
from ..score import score as scoreop
from ..io import wgt as wgtop
from ..validate.paper import auc as validate_paper_auc

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

    parser.add_argument('--model-boost',
                        help='')

    parser.add_argument('--model-boost-add',
                        help='')

    parser.add_argument('--methods',
                        help='')

    # use --para-regon
    # parser.add_argument('--regon-ss',
    #                    help='')

    parser.add_argument('--fpheno',
                        help='')

    parser.add_argument('--fcv',
                        help='')

    parser.add_argument('--cvn',
                        help='')

    parser.add_argument('--phes',
                        help='')

    parser.add_argument('--phes-female',
                        help='')

    parser.add_argument('--para-cand',
                        help='')

    parser.add_argument('--para-com',
                        help='')

    parser.add_argument('--para-regon',
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args

def main():


    args = argument()
    methods = args.methods.split(' ')
    phes = args.phes.split(' ')
    logger.debug('phes: {}'.format(phes))
    phes_female = args.phes_female.split(' ')
    if not set(phes_female).issubset(set(phes)):
        raise RuntimeError('Wrong phes_female', phes_female)

    boost_modeld = {
        'boosting_nonadd': args.model_boost,
        'boosting_add': args.model_boost_add
    }

    #paper(pathlib.Path(args.dout), pformats, pathlib.Path(args.dresult),
     #     pathlib.Path(args.ddata),
     #     methods, boost_modeld,
     #     pathlib.Path(args.fpheno), phes, phes_female,
     #     pathlib.Path(args.fcv), int(args.cvn),
     #     args.para_cand, args.para_regon, args.para_com,
     #     )


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
