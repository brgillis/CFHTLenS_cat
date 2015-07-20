#!/usr/bin/env python

import matplotlib
import os
import sys

import astropy.io.ascii as ascii
import icebrgpy.bins_funcs as bf
import matplotlib.pyplot as pyplot
import numpy as np
import subprocess as sbp


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


m_sun = 1.9891e30 # kg
pc = 149597870700. * 648000. / np.pi # m

def main(argv):
    
    # Magic values
    
    data_table_name = "/home/brg/git/CFHTLenS_cat/Data/alpha_cache_r_weighted.dat"
    
    paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    z_bins_min = 0.4
    z_bins_max = 1.5
    num_z_bins = 110
    z_bins_log = False
    
    z_bins_skip = 10
    
    indices_from_keys = {}
    keys_from_indices = {}
    col_i = 0
    
    col_names = []
    col_offsets = []
    col_factors = []
    
    indices_from_keys["z"] = col_i
    keys_from_indices[col_i] = "z"
    col_i += 1
    col_names.append("x_2")
    col_offsets.append(0)
    col_factors.append(1)
    
    indices_from_keys["mag"] = col_i
    keys_from_indices[col_i] = "mag"
    col_i += 1
    col_names.append("x_1")
    col_offsets.append(0.005)
    col_factors.append(1)
    
    indices_from_keys["alpha"] = col_i
    keys_from_indices[col_i] = "alpha"
    col_i += 1
    col_names.append("y")
    col_offsets.append(0)
    col_factors.append(1)
    
    num_cols = col_i

    # Load in the table
    try:
        data_table = ascii.read(data_table_name)
    except:
        print("ERROR: Table " + data_table_name + " cannot be read.")
        return
    
    cols = []
    binned_cols = []
           
    for col_i in xrange(num_cols):
        cols.append(np.array(data_table[col_names[col_i]])+col_offsets[col_i]*col_factors[col_i])
        binned_cols.append([])
        
    
    # Set up bins
    if(z_bins_log):
        z_bins, _z_bins_mids = bf.setup_log_bins(z_bins_min, z_bins_max, num_z_bins)
    else:
        z_bins, _z_bins_mids = bf.setup_lin_bins(z_bins_min, z_bins_max, num_z_bins)
    
    # Bin all values in the table
    for col_i in xrange(num_cols):
        for _ in xrange(num_z_bins):
            binned_cols[col_i].append([])

    for elems in zip(*cols):
        
        z = elems[indices_from_keys["z"]]
        
        z_i = bf.get_bin_index(z+0.001, z_bins)
        if(z_i<0): continue
        
        for col_i in xrange(num_cols):
            binned_cols[col_i][z_i].append(elems[col_i])
        
    # Convert all lists to arrays      
    for col_i in xrange(num_cols):
        for z_i in xrange(num_z_bins):
            binned_cols[col_i][z_i] = np.array(binned_cols[col_i][z_i])
        
    # Do the plotting now
    
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"$r'$ magnitude",labelpad=20)
    ax.set_xlim(20,24.7)
    ax.set_ylabel(r"$\alpha$",labelpad=20)
    ax.set_ylim(0,6)
    
    for z_i in xrange(0, num_z_bins, z_bins_skip):
        z = z_bins[z_i][0]
        
        label = "z="+str(z)
        
        ax.plot(binned_cols[indices_from_keys["mag"]][z_i],binned_cols[indices_from_keys["alpha"]][z_i],label=label)
        
    ax.legend(loc='upper right')
    
    
    # Save the figure
    outfile_name = os.path.splitext(data_table_name)[0] + ".eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)