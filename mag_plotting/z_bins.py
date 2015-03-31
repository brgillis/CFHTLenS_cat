'''
Created on 30 Mar 2015

@author: brg
'''

import numpy as np

def setup_z_bins(z_min,z_max,z_step):
        
    # Get the number of bins and buffer points to test
    num_z_bins = int(np.round((z_max-z_min)/z_step))
    
    # Get arrays of the mins, maxes, and mids for each bin
    z_mins = np.linspace(start=z_min,
                              stop=z_max,
                              num=num_z_bins,
                              endpoint=False)
    z_maxes = z_mins + z_step
    z_mids = (z_mins+z_maxes)/2.
    z_bins =  np.transpose(np.vstack((z_mins,z_maxes)))
    
    return num_z_bins, z_bins, z_mids
    
def get_bin_index(z, bins):
    ibest = np.argmax(z<bins[:,1])
    if(ibest>0):
        return ibest
    if((z>=bins[0,0]) and (z<bins[0,1])):
        return 0
    else:
        return -1