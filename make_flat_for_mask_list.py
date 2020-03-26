#!/usr/bin/env python


import numpy as np
import os


import sys
import parameters
import bfptc.filehandlers as filehandlers

if __name__ == "__main__" :
    pairs_filename = sys.argv[1]
    tuple_name = sys.argv[2]
    output_list_name = sys.argv[3]
    # the tuple contains time stamps, which are not easily converted into actual file names
    # so we cook up here the lookup table
    f = open(pairs_filename)
    fits_names = {}
    amps = None
    for i,l in enumerate(f.readlines()) :
        for fits in l.split() :
            im = filehandlers.SlacBot(fits, None)
            fits_names[im.time_stamp()] = fits
            if amps is None : amps = im.segment_ids()
    print("associated time stamps and fits file names")
    nt = np.load(tuple_name)
    nt = nt.view(np.recarray)
    nt = nt[(nt.i ==0) & (nt.j == 0)]
    chip = os.path.basename(os.path.dirname(os.path.abspath(tuple_name)))
    maxmu_amp = {amp:parameters.max_nonlin_value(chip, amp) for amp in amps}
    maxmu = np.array([maxmu_amp[x.ext] for x in nt])
    #nt = nt[(nt.mu1<maxmu) & (nt.mu1>0.5*maxmu)]
    nt = nt[(nt.mu1<maxmu)]
    time_stamps = set(nt.t1).union(set(nt.t2))
    output_list = open(output_list_name, "w")
    for time_stamp in time_stamps :
        output_list.write('%s\n'%fits_names[time_stamp])
















