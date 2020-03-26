#!/usr/bin/env python


#from ptc_utils import eval_nonlin
import numpy as np
import pylab as pl
import os

import scipy.interpolate as interp
import pickle



def plot_nonlin(nt, nonlin, chip, figname = None, liny = False, cut_high_tail=None) :
    amps = set(nt.ext)
    #f, all_axes = pl.subplots(4,4, sharex=True, sharey=True)
    f, all_axes = pl.subplots(4,4, figsize=(8,9), sharex = True)
    f.suptitle("Nonlinearity %s ref/CCD vs CCD"%chip)
    f.subplots_adjust(hspace=0.1, wspace=0.1)
    for k,amp in enumerate(amps):
        ntamp = nt[nt.ext ==  amp]
        ax = f.axes[k]
        # lin only
        # ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0)) 
        if not liny : ax.set_xscale('log')
        x = ntamp.mu1
        index = x.argsort()
        y = ntamp.d1[index]
        x = x[index]
        if cut_high_tail is not None :
            cut_value = x[-1]*cut_high_tail
            index = x<cut_value
            y = y[index]
            x = x[index]
        y = y * x.mean()/y.mean()
        points, = ax.plot(x, y/x-1, '.')
        if nonlin is not None :
            model = interp.splev(x, nonlin[amp])
            ax.plot(x, model/x, '-', color=points.get_color())
        if k == 12 :
            ax.set_xlabel("$\mu (ADU)$",fontsize='x-large')
            ax.set_ylabel("$diode/\mu -1 + Cst$",fontsize='x-large')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    #f.subplots_adjust(hspace=0, wspace=0)
    #pl.setp([a.get_xticklabels() for a in all_axes[1:4,1:4].flat], visible=False)
    #pl.tight_layout()
    if figname is not None : pl.savefig(figname)

# I don't know who wrote it...
def mad(data, axis=0, scale=1.4826):
    """
    Median of absolute deviation along a given axis.  

    Normalized to match the definition of the sigma for a gaussian
    distribution.
    """
    if data.ndim == 1:
        med = np.ma.median(data)
        ret = np.ma.median(np.abs(data-med))
    else:
        med = np.ma.median(data, axis=axis)
        if axis>0:
            sw = np.ma.swapaxes(data, 0, axis)
        else:
            sw = data
        ret = np.ma.median(np.abs(sw-med), axis=0)
    return scale * ret


def fit_nonlin_corr(xin, yclapin, knots = 20, loop = 20, verbose = False, fullOutput=False):
    """
    xin : the data to be "linearized"
    yclapin : the (hopefully) linear reference
    returns a  spline that can be used for correction uisng "scikit.splev"
    if full_output==True, returns spline,x,y  (presumably to plot)
    """
    # do we need outlier rejection ?
    # the xin has to be sorted, although the doc does not say it....
    index = xin. argsort()
    x = xin[index]
    yclap = yclapin[index]
    chi2_mask = np.isfinite(yclap) # yclap = nan kills the whole thing
    xx = x
    yyclap = yclap
    for i in range(loop):
        xx = xx[chi2_mask]
        # first fit the scaled difference between the two channels we are comparing
        yyclap = yyclap[chi2_mask]
        length = xx[-1]-xx[0]
        t = np.linspace(xx[0]+1e-5*length, xx[-1]-1e-5*length, knots)
        s = interp.splrep(xx, yyclap, task=-1, t=t)        
        model = interp.splev(xx, s)     # model values
        res = model - yyclap
        sig = mad(res)
        res = np.abs(res)
        if (res> (5 * sig)).sum()>0 : # remove one at a time
            chi2_mask = np.ones(len(xx)).astype(bool)
            chi2_mask[np.argmax(res)] = False
            continue
        else : break
    # enforce the fit to got through (0,0) and make sure that 0 is inside
    # the definition domain of the spline.
    # print('means yy ',yyclap.mean(),' xx ', xx.mean())
    # print ('ymod[0]/xmod[0] ', interp.splev(xx[0],s)/xx[0],xx[0])
    old_der = yyclap.mean()/xx.mean()
    nx = len(xx)
    fact = 1
    nadd = nx//2
    fit_val = interp.splev(xx[0],s)
    fake_x = np.linspace(-xx[0]*fact, xx[0]*fact, nadd)
    fake_y = np.linspace(-fit_val*fact, fit_val*fact, nadd)
    xx = np.hstack((fake_x , xx))
    yyclap = np.hstack((fake_y , yyclap))
    t = np.linspace(xx[0]+1e-5*length, xx[-1]-1e-5*length, knots)
    s = interp.splrep(xx, yyclap, task=-1, t=t[1:-2])
    # normalize to "no change" at x->0
    der0 = interp.splev(0., s, 1)
    norm_fact = 1./der0
    # print("n before/after %d/%d"%(nx,len(xx)))
    norm_fact = yyclap.mean()/xx.mean()
    # print('comparison old_fact ', old_der, ' new_fact ',norm_fact)
    yyclap_norm = yyclap / norm_fact
    # model only the residual to the identity
    s = interp.splrep(xx, yyclap_norm - xx , task=-1, t=t)

    model = interp.splev(xx, s) + xx    # model values
    # compute gain residuals
    mask = (yyclap_norm != 0)
    print("der0 = %f, val0 = %f"%(1+interp.splev(0., interp.splder(s)),interp.splev(0.,s)),"nonlin gain residuals : %g"%(model[mask]/yyclap_norm[mask]-1).std())
    if verbose :     
        print("fit_nonlin loops=%d sig=%f res.max = %f"%(i,sig, res.max()))
    if fullOutput :
        return s, xx, yyclap_norm
    return s

def eval_nonlin(tuple, knots = 20, verbose = False, fullOutput=False, ref_name='c'):
    """
    it will be faster if the tuple only contains the variances
    return value: a dictionnary of correction spline functions (one per amp)
    """
    amps = np.unique(tuple['ext'].astype(int))
    res={}
    if fullOutput:
        x = {}
        y = {}
    rname = ref_name+'1'
    for i in amps :
        t = tuple[tuple['ext'] == i]
        clap = np.hstack((t[ref_name+'1'],t[ref_name+'2']))
        mu = np.hstack((t['mu1'],t['mu2']))
        if fullOutput :
            res[i], x[i], y[i] = fit_nonlin_corr(mu,clap, knots=knots, verbose=verbose, fullOutput=fullOutput)
        else :
            res[i] = fit_nonlin_corr(mu,clap, knots=knots, verbose=verbose, fullOutput=fullOutput)
    if fullOutput:
        return res,x,y
    else :
        return res






import sys
import parameters
import pickle


if __name__ == "__main__" :
    tuple_name = sys.argv[1]
    nt = np.load(tuple_name)
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple_name)))
    nt = nt.view(np.recarray)
    # the filter2.pkl should sit in the same directory as this script
    where = os.path.dirname(sys.argv[0])
    filter2 = pickle.load(open(where+'/filter2.pkl','rb'))
    nt = nt[(nt.i ==0) & (nt.j == 0)]
    # cut on value of filter2 (discussion on #eochar on that)
    f2 = [filter2[x]=='ND_OD0.01' for x in nt.t1]
    nt_proc = nt[f2]
    # loop on channels because the max value cut depends on channel
    amps = set(nt_proc.ext)
    maxmu = {amp:parameters.max_nonlin_value(chip, amp) for amp in amps}
    maxmu_array =  [maxmu[i] for i in nt_proc.ext]
    nt_cut_mu = nt_proc[nt_proc.mu1 < maxmu_array]
    res = eval_nonlin(nt_cut_mu, knots = parameters.spline_knots_for_nonlin, ref_name='d')
    pickle.dump(res,open(sys.argv[2],"wb"))
    pl.ioff()
    plot_nonlin(nt_cut_mu, res, chip, 'nonlin.png')
    plot_nonlin(nt, None, chip, 'nonlin_full.png', liny = True, cut_high_tail=0.9)
















