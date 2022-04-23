    
from altaipony.lcio import from_mast, from_path, FlareLightCurve
import numpy as np
import time
import copy
import pandas as pd

from altaipony.altai import find_iterative_median
from altaipony.utils import sigma_clip

from scipy import constants as c

def blackbodyw(wavelength, temperature, response): 
    pre_factor = 2* (c.h * c.c**2)/ wavelength**5 
    exponential_factor = 1. /(np.exp((c.h*c.c) / (wavelength * c.k * temperature))-1)
    return pre_factor * exponential_factor * response #fake convolution with response

from scipy import integrate
#energy integration
def energy_filter(wavelength, T ,response):
    y1 = blackbodyw(wavelength, T, response) 
    return integrate.trapz(y1)

def energy_filter2(response):
    y1 = response 
    return integrate.trapz(y1)

def lcbin(time, flux, binwidth=0.5, nmin=4, time0=None,
        robust=False, tmid=False):
    '''
    Calculate average flux and error in time bins of equal width.
    The default bin width is equivalent to one CHEOPS orbit in units of days.
    To avoid binning data on either side of the gaps in the light curve due to
    the CHEOPS orbit, the algorithm searches for the largest gap in the data
    shorter than binwidth and places the bin edges so that they fall at the
    centre of this gap. This behaviour can be avoided by setting a value for
    the parameter time0.
    The time values for the output bins can be either the average time value
    of the input points or, if tmid is True, the centre of the time bin.
    If robust is True, the output bin values are the median of the flux values
    of the bin and the standard error is estimated from their mean absolute
    deviation. Otherwise, the mean and standard deviation are used.
    The output values are as follows.
    * t_bin - average time of binned data points or centre of time bin.
    * f_bin - mean or median of the input flux values.
    * e_bin - standard error of flux points in the bin.
    * n_bin - number of flux points in the bin.
    :param time: time
    :param flux: flux (or other quantity to be time-binned)
    :param binwidth:  bin width in the same units as time
    :param nmin: minimum number of points for output bins
    :param time0: time value at the lower edge of one bin
    :param robust: use median and robust estimate of standard deviation
    :param tmid: return centre of time bins instead of mean time value
    :returns: t_bin, f_bin, e_bin, n_bin
    '''
    if time0 is None:
        tgap = (time[1:]+time[:-1])/2
        gap = time[1:]-time[:-1]
        j = gap < binwidth
        gap = gap[j]
        tgap = tgap[j]
        time0 = tgap[np.argmax(gap)]
        time0 = time0 - binwidth*np.ceil((time0-min(time))/binwidth)
        
    n = np.int(1+np.ceil(np.ptp(time)/binwidth))
    r = (time0,time0+n*binwidth)
    
    n_in_bin,bin_edges = np.histogram(time,bins=n,range=r)
    bin_indices = np.digitize(time,bin_edges)
    
    t_bin = np.zeros(n)
    f_bin = np.zeros(n)
    e_bin = np.zeros(n)
    n_bin = np.zeros(n, dtype=np.int)
    
    for i,n in enumerate(n_in_bin):
        if n >= nmin:
            j = bin_indices == i+1
            n_bin[i] = n
            if tmid:
                t_bin[i] = (bin_edges[i]+bin_edges[i+1])/2
            else:
                t_bin[i] = np.mean(time[j])
            if robust:
                f_bin[i] = np.median(flux[j])
                e_bin[i] = 1.25*np.mean(abs(flux[j] - f_bin[i]))/np.sqrt(n)
            else:
                f_bin[i] = np.mean(flux[j])
                e_bin[i] = np.std(flux[j])/np.sqrt(n-1)
    j = (n_bin >= nmin)
    return t_bin[j], f_bin[j], e_bin[j], n_bin[j]
    


def estimate_detrended_noise(flc, mask_pos_outliers_sigma=2.5, std_window=100, padleft=3, padright=10):

    flcc = copy.deepcopy(flc)
    flcc = flcc.find_gaps()

    for (le, ri) in flcc.gaps:

        flcd = copy.deepcopy(flcc[le:ri])
        mask = sigma_clip(flcd.detrended_flux, max_sigma=mask_pos_outliers_sigma)#, longdecay=2)

        flcd.detrended_flux[~mask] = np.nan
        # apply rolling window std and interpolate the masked values
        flcd.detrended_flux_err[:] = pd.Series(flcd.detrended_flux).rolling(std_window,
                                                                 center=True,
                                                                 min_periods=1).std().interpolate()
         
        # and refine it:
        flcd = find_iterative_median(flcd)
        
        # make a copy first
        filtered = copy.deepcopy(flcd.detrended_flux)
        
        # get right bound of flux array
        tf = filtered.shape[0]

        # pick outliers
        mask = sigma_clip(filtered, max_sigma=mask_pos_outliers_sigma)#, longdecay=2)


        # apply rolling window std and interpolate the masked values
        flcc.detrended_flux_err[le:ri]= pd.Series(filtered).rolling(std_window,
                                                                 center=True,
                                                                 min_periods=1).std().interpolate()
    return flcc
    


