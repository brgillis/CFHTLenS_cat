#!/usr/bin/env python

import os
import subprocess as sbp
import sys

import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import bins_funcs as bf
import columns_reader as cr

m_sun = 1.9891e30 # kg
pc = 149597870700. * 648000. / np.pi # m

sqdeg_per_sqrad = (180/np.pi)**2

def main(argv):
    
    # Magic values
    
    data_table_name = "/home/brg/git/CFHTLenS_cat/Data/ex_count_cache_r_weighted.dat"
    
    paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    z_bins_min = 0.4
    z_bins_max = 1.5
    num_z_bins = 110
    z_bins_log = False
    
    z_bins_skip = 10
    
    reader = cr.columns_reader()
    
    reader.add("z","x_2",offset=0.001)
    reader.add("mag","x_1",offset=0.005)
    reader.add("count","y")
        
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
    scale = 1/(sqdeg_per_sqrad)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"$r'$" + " magnitude",labelpad=20)
    ax.set_xlim(20,24.7)
    
    ax.set_ylabel(r"Count per deg$^2$ per unit magnitude",labelpad=20)
    ax.set_yscale("log", nonposy='clip')
    ax.set_ylim(1e4*scale,1.5e8*scale)
    
    plotted_y_cols = []

    cols = reader.read(data_table_name)
    
    binned_cols = []
           
    for col_i in xrange(reader.num_cols()):
        binned_cols.append([])
        
    
    # Set up bins
    if(z_bins_log):
        z_bins, _z_bins_mids = bf.setup_log_bins(z_bins_min, z_bins_max, num_z_bins)
    else:
        z_bins, _z_bins_mids = bf.setup_lin_bins(z_bins_min, z_bins_max, num_z_bins)
    
    # Bin all values in the table
    for col_i in xrange(reader.num_cols()):
        for _ in xrange(num_z_bins):
            binned_cols[col_i].append([])

    for elems in zip(*cols):
        
        z = elems[reader.index("z")]
        
        z_i = bf.get_bin_index(z+0.001, z_bins)
        if(z_i<0): continue
        
        for col_i in xrange(reader.num_cols()):
            binned_cols[col_i][z_i].append(elems[col_i])
        
    # Convert all lists to arrays      
    for col_i in xrange(reader.num_cols()):
        for z_i in xrange(num_z_bins):
            binned_cols[col_i][z_i] = np.array(binned_cols[col_i][z_i])
        
    # Do the plotting now
    
    
    for z_i in xrange(0, num_z_bins, z_bins_skip):
        z = z_bins[z_i][0]
    
        if("unweighted" in data_table_name):
            lstyle = '--'
            label = None
        else:
            lstyle = '-'
            label = r"z$\geq$"+str(z)
        
        ax.plot(binned_cols[reader.index("mag")][z_i],binned_cols[reader.index("count")][z_i]*scale,
                lstyle,label=label)
        
        plotted_y_cols.append(binned_cols[reader.index("count")][z_i]*scale)
            
    ax.legend(loc='lower right')
        
        
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