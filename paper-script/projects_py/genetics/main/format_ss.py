"""
For assoc
"""


import argparse
from logging import getLogger
from ...system.logger import logger_setting
# import decimal
from decimal import Decimal

#from ..io_genot import snv as snvio
from ..snv import snvio as snvio


logger = getLogger(__name__)


def argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fout',
                        default=None,
                        help='output file')

    parser.add_argument('--add_a2',
                        action='store_true',
                        help='')

    parser.add_argument('--add_col_p',
                        action='store_true',
                        help='')

    parser.add_argument('--fassoc',
                        default=None,
                        help='')

    args = parser.parse_args()

    logger.info('args {}'.format(args))

    return args


def add_col_p(fout, fassoc):
    def str_from_log10_decimal(x):
        # https://stackoverflow.com/questions/36092203/why-is-1e400-not-an-int
        dec = Decimal(10)**(Decimal(-1) * Decimal(x))
        if dec < 10**(-4):
            num_str = '{:.6e}'.format(dec)
        else:
            num_str = '{:.6f}'.format(dec)
        return num_str
        # return str(decimal.Decimal(x).log10())
        # ng
        # python float cannot handle 1e-350 as well
        # return math.log10(float(x))

    assoc = snvio.load_assoc(fassoc)
    assoc.loc[:, 'P'] = assoc['LOG10_P'].apply(lambda x: str_from_log10_decimal(x))

    assoc.to_csv(fout, sep='\t', index=False)


def add_a2(fout, fassoc):
    def get_a2(x):
        if x['A1'] == x['ALT']:
            return x['REF']
        elif x['A1'] == x['REF']:
            return x['ALT']
        else:
            raise RuntimeError('A1 if not either ALT or REF')

    assoc = snvio.load_assoc(fassoc)

    # remove #
    assoc.columns = assoc.columns.str.replace('#', '')

    if (assoc['ALT']!=assoc['A1']).any():
        raise RuntimeError('Use "omit-ref" option in plink, or comment out this assertion')

    logger.debug('assoc\n{}'.format(assoc.head()))
    assoc.loc[:, 'A2'] = assoc[['REF', 'ALT', 'A1']].apply(lambda x: get_a2(x), axis='columns')

    assoc.to_csv(fout, sep='\t', index=False)


def main():
    args = argument()

    if args.add_a2:
        # logger.info('join_snv_hm3]')
        add_a2(args.fout, args.fassoc)

    if args.add_col_p:
        # logger.info('join_snv_hm3]')
        add_col_p(args.fout, args.fassoc)


if __name__ == '__main__':
    logger_setting.setting()
    main()
    logger.debug('Done!')
