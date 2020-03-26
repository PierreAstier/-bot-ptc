#!/usr/bin/env python

import numpy as np
import os.path
import matplotlib.pyplot as pl
import sys


def plot_npix(fullnt, chip, mu_max, figname=None) :
    """
    plot npix%mu
    """
    f, all_axes = pl.subplots(4,4, figsize=(10,10))
    f.subplots_adjust(hspace=0.1, wspace=0.1)
    f.suptitle('%s npix vs mu '%chip)
    nt = fullnt[(fullnt.i == 0) & (fullnt.j == 0)]
    channels = set(nt.ext)
    for k,channel in enumerate(channels) :
        nte = nt[(nt.ext == channel) & (nt.mu1<mu_max)]
        if k == 0 : channel0 = channel
        ax = f.axes[channel-channel0]
        x = (nte.mu1+nte.mu2)*0.5
        y = nte.npix
        ax.plot(x,y,'.')
        if k%4 == 0 : ax.set_ylabel('npix', fontsize='x-large')
        if k//4 == 3 : ax.set_xlabel('$\mu$', fontsize='x-large')
    if figname is not None : pl.savefig(figname)


if __name__ == "__main__" :
    tuple = sys.argv[1]
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple)))
    nt = np.load(tuple)
    nt = nt.view(np.recarray)
    plot_npix(nt, chip, mu_max=1.2e5, figname=None)
    pl.show(block=True)
