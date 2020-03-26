#!/usr/bin/env python

import pickle
import numpy as np
import ptc_utils
import os.path
import bfptc.ptcfit as ptcfit
import matplotlib.pyplot as pl
import sys
import os

def plot_ptc(fits, chip, i=0,j=0, figname = None, plot_model = True) :
    #fig = pl.figure(figsize=(8,8))
    #f, all_axes = pl.subplots(4,4, sharex=True, sharey=True)
    f, all_axes = pl.subplots(4,4, figsize=(8,8))
    f.suptitle('%s $C_{%d%d}/\mu$'%(chip,i,j))
    for channel,fit in fits.items():
        ax = f.axes[channel-1]
        index = (fit.sqrt_w[:,i,j] != 0)
        mu = fit.mu[index]
        ax.plot(mu, fit.cov[index,i,j]/mu, '.')
        if plot_model :
            model = fit.eval_cov_model()[index,i,j]
            ax.plot(mu, model/mu,'-')
    #f.subplots_adjust(hspace=0, wspace=0)
    #pl.setp([a.get_xticklabels() for a in all_axes[1:4,1:4].flat], visible=False)
    if figname is not None : pl.savefig(figname)

def plot_ptc_slice(nt, chip, i, j, figname = None) :
    """
    averages covariances over some range in i and j (which should be slices)
    the mu cuts should have been applied upstream
    
    """
    #f, all_axes = pl.subplots(4,4, sharex=True, sharey=True)
    f, all_axes = pl.subplots(4,4, figsize=(8,8))
    f.suptitle('%s average $C_{ij}/\mu$ for i,j in [%d,%d[ [%d,%d[ '%
               (chip,i.start,i.stop,j.start,j.stop))
    exts = set(nt.ext)
    for k,ext in enumerate(exts) :
        ax = f.axes[k]
        nte = nt[nt.ext == ext]
        cov,vcov,mu = ptcfit.make_cov_array(nte,r=max(i.stop,j.stop))
        y = cov[:,i,j].mean(axis=(1,2))
        ax.plot(mu, y/mu, '.')
    #f.subplots_adjust(hspace=0, wspace=0)
    #pl.setp([a.get_xticklabels() for a in all_axes[1:4,1:4].flat], visible=False)
    if figname is not None : pl.savefig(figname)

def plot_a_vs_dist(fits, chip, figname = None) :
    f, all_axes = pl.subplots(4,4, figsize=(10,10),sharex=True, sharey=True)
    f.subplots_adjust(hspace=0, wspace=0)
    f.suptitle('%s : $a_{ij}$  vs $\sqrt{i^2+j^2}$'%chip)
    for channel,fit in fits.items():
        ax = f.axes[channel-1]
        a = fit.get_a()
        i,j=np.indices(a.shape)
        upper = (i>=j).ravel()
        r = np.sqrt(i**2+j**2).ravel()
        a = a.ravel()
        ax.plot(r[upper],a[upper], marker='o', color='b', linestyle='none', label='$i>=j$')
        ax.plot(r[~upper],a[~upper], marker='o', color='r', linestyle='none', label='$i< j$')
        # negative ones
        ax.plot(r[upper],-a[upper], marker='v', color='b', linestyle='none', label='$i>=j$ $a<0$')
        ax.plot(r[~upper],-a[~upper], marker='v', color='r', linestyle='none', label='$i< j$ $a<0$')
        ax.set_yscale('log')
        if channel == 1 : ax.legend(loc='best')
    #pl.tight_layout()
    if figname is not None : pl.savefig(figname)

def plot_cov_average(nt, chip, mu_min, mu_max, maxr=20, figname=None) :
    """
    averages covariances in a signal level range and plots 
    the result as an image
    """
    f, all_axes = pl.subplots(4,4, figsize=(11,10))
    f.subplots_adjust(hspace=0.1, wspace=0.1)
    f.suptitle('%s average covariances for %f < $\mu$ < %f '%(chip,mu_min,mu_max))
    channels = set(nt.ext)
    for channel in channels :
        nte = nt[nt.ext == channel]
        ax = f.axes[channel-1]
        cov,vcov,mu = ptcfit.make_cov_array(nte,r=maxr)
        index_mu = (mu>mu_min) & (mu<mu_max)
        cov_mean = cov[index_mu,:,:].mean(axis=0)
        if channel == 1 :
            m = np.median(cov_mean)
            s = np.median(np.abs(m-cov_mean))
            vmin = m-5*s
            vmax = m+5*s
        im = ax.imshow(cov_mean.T,origin='lower',vmin=vmin,vmax=vmax)
        ax.set_xlabel('i', fontsize='x-large')
        ax.set_ylabel('j', fontsize='x-large')
    # draw the colorbar on the RHS
    f.subplots_adjust(right=0.9)
    cbar_ax = f.add_axes([0.91, 0.1, 0.03, 0.8])
    f.colorbar(im,cax = cbar_ax)
    if figname is not None : pl.savefig(figname)

def plot_noise(fits, chip, maxr=20, figname=None) :
    """
    plot fitted noise image
    """
    f, all_axes = pl.subplots(4,4, figsize=(11,10))
    f.subplots_adjust(hspace=0.1, wspace=0.1)
    f.suptitle('%s Noise covariances (el$^2$)'%(chip))
    for channel,fit in fits.items():
        ax = f.axes[channel-1]
        n = fit.get_noise()
        if channel == 1 :
            m = np.median(n)
            s = np.median(np.abs(m-n))
            vmin = m-20*s
            vmax = m+20*s
        im = ax.imshow(n.T,origin='lower',vmin=vmin,vmax=vmax)
        ax.set_xlabel('i', fontsize='x-large')
        ax.set_ylabel('j', fontsize='x-large')
    # draw the colorbar on the RHS
    f.subplots_adjust(right=0.9)
    cbar_ax = f.add_axes([0.91, 0.1, 0.03, 0.8])
    f.colorbar(im,cax = cbar_ax)
    if figname is not None : pl.savefig(figname)

    

if __name__ == "__main__" :
    tuple = sys.argv[1]
    output = sys.argv[2]
    print("fitting %s"%tuple)
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple)))
    nt = np.load(tuple)
    nt = nt.view(np.recarray)
    params = ptc_utils.load_params()
    params.subtract_distant_value = False # Not sure
    params.maxmu_el = 9e4
    params.maxmu = 8e4
    lfit = ptc_utils.load_data(nt,params)
    pl.ioff()
    # plots without fits
    mu_min = 0.5*params.maxmu
    mu_max = params.maxmu
    plot_noise(lfit, chip, maxr=20, figname='noise_cov.png')
    plot_cov_average(nt, chip, mu_min, mu_max, maxr=20, figname='all_cov.png')
    plot_ptc_slice(nt[(nt.mu1<params.maxmu)], chip, slice(10,20), slice(0,1), figname ='cov_serial10.png')
    plot_ptc_slice(nt[(nt.mu1<params.maxmu)], chip, slice(20,30), slice(0,1), figname ='cov_serial20.png')
    plot_ptc(lfit,chip, 0,0,'ptc_data.png',plot_model = False)
    plot_ptc(lfit,chip, 1,0,'C10_data.png',plot_model = False)
    plot_ptc(lfit,chip, 0,1,'C01_data.png',plot_model = False)
    # now fit:
    for k,c in lfit.items():
        if c is None : continue
        print('fitting channel %d'%k)
        c.fit()
    f = open(output,'wb')
    pickle.dump(lfit,f)
    f.close()
    plot_noise(lfit, chip, maxr=20, figname='noise_cov.png')
    plot_ptc(lfit, chip, 0,0, 'ptcfit.png')
    plot_ptc(lfit, chip, 1,0, 'C10fit.png')
    plot_ptc(lfit, chip, 0,1, 'C01fit.png')
    plot_a_vs_dist(lfit, chip, figname = 'a_vs_dist.png')
    la = []
    for _,c in lfit.items():
        la.append(c.get_a())
    la = np.array(la)
    am = np.median(la,axis = 0)
    maxr = min(am.shape[0],3)
    print('median  a:\n',am[:maxr,:maxr].T)
    print('a00 : ',la[:,0,0])
