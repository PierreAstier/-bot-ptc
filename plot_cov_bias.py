#!/usr/bin/env python

import numpy as np
import ptcfit
import os.path
import matplotlib.pyplot as pl
import sys
import os

def plot_cov_average(nt, chip, maxr=30, figname=None) :
    """
    averages covariances in a signal level range and plots 
    the result as an image
    """
    f, all_axes = pl.subplots(4,4, figsize=(11,10))
    f.subplots_adjust(hspace=0.1, wspace=0.1)
    f.suptitle('%s average covariances for biases'%(chip))
    channels = set(nt.ext)
    for k,channel in enumerate(channels) :
        nte = nt[nt.ext == channel]
        ax = f.axes[channel-1]
        cov,vcov,mu = ptcfit.make_cov_array(nte,r=maxr)
        cov_mean = cov.mean(axis=0)
        if channel == 1 :
            m = np.median(cov_mean)
            s = np.median(np.abs(m-cov_mean))
            vmin = m-5*s
            vmax = m+5*s
        im = ax.imshow(cov_mean.T,origin='lower',vmin=vmin,vmax=vmax)
        if k//4 == 3 : ax.set_xlabel('i', fontsize='x-large')
        if k%4 == 0 : ax.set_ylabel('j', fontsize='x-large')
    # draw the colorbar on the RHS
    f.subplots_adjust(right=0.9)
    cbar_ax = f.add_axes([0.91, 0.1, 0.03, 0.8])
    f.colorbar(im,cax = cbar_ax)
    if figname is not None : pl.savefig(figname)



if __name__ == "__main__" :
    tuple = sys.argv[1]
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple)))
    nt = np.load(tuple)
    nt = nt.view(np.recarray)
    pl.ioff()
    plot_cov_average(nt, chip, maxr=30, figname='cov_bias.png')

