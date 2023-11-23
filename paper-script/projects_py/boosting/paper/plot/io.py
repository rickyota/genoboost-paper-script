

from logging import getLogger

logger = getLogger(__name__)


def savefig_png_svg(fig, fname, dpi=None):
    fig.savefig(fname.parent / (fname.name + '.png'))
    # fig.savefig(fname.parent / (fname.name + '.pdf'))

    if dpi is None:
        fig.savefig(fname.parent / (fname.name + '.svg'))
    else:
        fig.savefig(fname.parent / (fname.name + '.svg'), dpi=dpi)
