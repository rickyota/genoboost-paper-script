
# print wgt path with the best para

import argparse
import pathlib
# from collections import defaultdict
from logging import getLogger

from ...system.logger import logger_setting
# from ...genetics.sample import sampleio as sampleio
from ..io import filename as iof
from . import paper

logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--ddata',
                        help='')

    # ./result/nonaddpgs/dataset.12/
    parser.add_argument('--dresult',
                        help='output dir')

    parser.add_argument('--method',
                        help='')

    parser.add_argument('--mid-path',
                        nargs="+",
                        help='')

    parser.add_argument('--phe',
                        help='')

    parser.add_argument('--cvi',
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


def print_fwgt_best(dresult, methodm, mid_path, phe, cvi, fymld):

    max_on = 'nagelkerke'
    withcov_at_maxon = True

    # load_manhattan_ssd(dresult, methodm, method_modeld, phe, cvi, fymld, max_on, withcov_at_maxon)

    regond = paper.load_regon(fymld)

    method_modeld = {methodm: mid_path}

    reg_on, va_on = paper.regon_vaon(methodm, regond, withcov_at_maxon)
    paras_best_cvs = paper.load_best_para_cvs(dresult, methodm, phe, withcov_at_maxon, reg_on, max_on,
                                              va_on, fymld['cand'], method_modeld, None, None, cvis=[cvi])

    # wgt = load_best_wgt(dresult, methodm, phe, cvi, paras_best_cvs, method_modeld)

    method = paper.methodm_to_method(methodm)
    # None if not boosting
    # mid_path = method_modeld[methodm]

    paras_best = paras_best_cvs[cvi]
    logger.debug('paras_best method phe: {}, {}, {}'.format(methodm, phe, paras_best))
    # should be methodm='boosting_nonadd'
    fwgt = iof.file_wgt(dresult, method, phe, cvi, paras_best, mid_path)

    nsnv = paras_best['nsnv_use']

    # pass to bash
    print(fwgt, nsnv)


def get_fwgt_best(ddata, dresult, methodm, mid_path,
                  phe, cvi,
                  para_cand, para_regon, para_com,
                  ):

    fymld = paper.get_fymld(ddata, para_cand, para_regon, para_com)

    print_fwgt_best(dresult, methodm, mid_path, phe, cvi, fymld,)


def main():

    args = argument()

    # method-> mid_path
    # ex. (kind), (kind, model), (kind, simu_type, model)
    # method_modeld = paper.arg_method_mid_path(args.method_model)
    # logger.debug('method_modeld {}'.format(method_modeld))

    # TODO: implement integrate

    get_fwgt_best(pathlib.Path(args.ddata), pathlib.Path(args.dresult),
                  args.method, tuple(args.mid_path),
                  args.phe, int(args.cvi),
                  args.para_cand, args.para_regon, args.para_com)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
