'''
Created on 30 Mar 2015

@author: brg
'''

import numpy as np

def setup_lin_bins(bins_min,bins_max,num_bins):
    
    step = float(bins_max-bins_min)/num_bins
    
    # Get arrays of the mins, maxes, and mids for each bin
    mins = np.linspace(start=bins_min,
                              stop=bins_max,
                              num=num_bins,
                              endpoint=False)
    maxes = mins + step
    mids = (mins+maxes)/2.
    bins =  np.transpose(np.vstack((mins,maxes)))
    
    return bins, mids

def setup_log_bins(bins_min,bins_max,num_bins):
    
    # Get arrays of the mins, maxes, and mids for each bin
    loglimits = np.linspace(start=np.log(bins_min),
                              stop=np.log(bins_max),
                              num=num_bins+1,
                              endpoint=True)
    limits = np.exp(loglimits)
    
    mins = limits[:-1]
    maxes = limits[1:]
    mids = (mins+maxes)/2.
    bins =  np.transpose(np.vstack((mins,maxes)))
    
    return bins, mids
    
def get_bin_index(x, bins):
    ibest = np.argmax(x<bins[:,1])
    if(ibest>0):
        return ibest
    if((x>=bins[0,0]) and (x<bins[0,1])):
        return 0
    else:
        return -1