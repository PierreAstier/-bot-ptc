#!/usr/bin/env python



import numpy as np
import parameters
import scipy.interpolate as interp

def apply_cuts(nt, mumax, pixmax):
    cuts = (nt.mu<mumax) & (nt.im<pixmax)
    return nt[cuts]


def fit_over(nt):
    amps = set(nt.amp.astype(int))
    res = {}
    for amp in amps:
        print('fitting amp %d'%amp)
        nta = nt[(nt.amp == amp)]
        x = nta.im
        index = np.argsort(x)
        x = x[index]
        y = nta.o1[index]
        nknots = parameters.nknots_for_cti_model
        t = np.linspace(x.min(), x.max(), nknots+2)
        s = interp.splrep(x,y,task=-1, t=t[1:-2])
        res[amp] = s
    return res

import matplotlib.pyplot as pl
import os

def plot_over(nt, models, title, figname = None) :
    amps = set(nt.amp)
    f, all_axes = pl.subplots(4,4,figsize=(9,9))
    f.suptitle(chip)
    for k,amp in enumerate(amps):
        print('plotting %d'%amp)
        nta = nt[(nt.amp == amp)]
        ax = f.axes[k]
        x = nta.im
        npoints = len(x)
        sampling = int(npoints/10000)
        ax.plot(x[::sampling], nta.o1[::sampling], ',')
        xc = np.linspace(x.min(), x.max(),100)
        ax.plot(xc, interp.splev(xc,models[amp]), 'k-')
        #ax.set_title('Amp %d'%amp)
    #print('calling tight_layout')
    #pl.tight_layout() 
    #print('saving figure')
    if figname is not None : pl.savefig(figname)
    

    
import sys
import pickle
if __name__ == "__main__" :
    tuple_name = sys.argv[1]
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple_name)))
    vendor = parameters.ccd_vendor(chip)
    output_name = sys.argv[2]
    nt = np.load(tuple_name)
    nt = nt.view(np.recarray)
    ntc = apply_cuts(nt, parameters.max_mu_for_over_fit, parameters.max_pix_for_over_fit)
    models = fit_over(ntc)
    pickle.dump(models, open(output_name,'wb'))
    pl.ioff()
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple_name)))
    plot_over(ntc, models, chip, "cti.png")
