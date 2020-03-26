spline_knots_for_nonlin = 10
nknots_for_cti_model = 10
max_mu_for_over_fit = 1.1e5
max_pix_for_over_fit = 1.1e5
degree_for_cti_model = 3
mumax_for_nonlin_fit = 1.1e5
ccd_vendor_dict={'R11':'e2v',
 'R12':'e2v',
 'R21':'e2v',
 'R22':'e2v',
 'R30':'e2v',
 'R00':'ITL',
 'R01':'ITL',
 'R02':'ITL',
 'R04':'ITL',
 'R10':'ITL',
 'R20':'ITL'}

def ccd_vendor(chip) :
    return ccd_vendor_dict[chip[:3]]

def max_nonlin_value(chip,channel) :
    if ccd_vendor(chip) == 'e2v' : return mumax_for_nonlin_fit
    # it is an ITL, lower for channel 9
    if channel == 9: return 0.8*mumax_for_nonlin_fit
    return  mumax_for_nonlin_fit
