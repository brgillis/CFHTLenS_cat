#!/usr/bin/env python

import os
import subprocess as sbp
import sys

import numpy as np
import astropy.io.ascii as ascii

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import bins_funcs as bf

def main(argv):
    
    # Magic values
    default_redshift_table_name = "/home/brg/git/CFHTLenS_cat/Data/spec_matched_tables/all_matched_galaxies.dat"
    default_z_spec_name = "Z"
    default_z_phot_name = "Z_B"
    default_z_bins_min = 0.2
    default_z_bins_max = 1.3
    default_num_z_bins = 11
    default_z_buffer_min = 0
    default_z_buffer_max = 0.3
    default_z_buffer_step = 0.05
    
    default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    # Load in values from command line if provided, otherwise use defaults
    
    num_args = len(argv)
    cur_arg = 2
    
    if(num_args) < cur_arg:
        redshift_table_name = default_redshift_table_name
    else:
        redshift_table_name = argv[cur_arg]
    
    # Load in the redshift table
    try:
        redshift_table = ascii.read(redshift_table_name)
    except:
        print("ERROR: Redshift table " + redshift_table_name + " cannot be read.")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_spec_name = default_z_spec_name
    else:
        z_spec_name = argv[cur_arg]
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_phot_name = default_z_phot_name
    else:
        z_phot_name = argv[cur_arg]
        
    # Try to load in the columns for spec and phot redshifts
    try:
        spec_zs = redshift_table[z_spec_name]
        phot_zs = redshift_table[z_phot_name]
        spec_and_phot_zs = np.transpose(np.vstack((spec_zs,phot_zs)))
    except:
        print("ERROR: Cannot read columns " + z_spec_name + " and " + z_phot_name + " from table " + redshift_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_bins_min = default_z_bins_min
    else:
        z_bins_min = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_bins_max = default_z_bins_max
    else:
        z_bins_max = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        num_z_bins = default_num_z_bins
    else:
        num_z_bins = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_buffer_min = default_z_buffer_min
    else:
        z_buffer_min = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_buffer_max = default_z_buffer_max
    else:
        z_buffer_max = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_buffer_step = default_z_buffer_step
    else:
        z_buffer_step = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        paper_location = default_paper_location
    else:
        paper_location = argv[cur_arg]
        
    # Get the number buffer points
    num_z_buffers = int(np.round((z_buffer_max-z_buffer_min)/z_buffer_step))+1
    
    # Set up the z bins arrays and data
    z_bins, z_bins_mids = bf.setup_lin_bins(bins_min=z_bins_min, bins_max=z_bins_max, num_bins=num_z_bins)
    
    # Get an array of buffer points
    z_buffers = np.linspace(start=z_buffer_min,
                              stop=z_buffer_max,
                              num=num_z_buffers,
                              endpoint=True)
    # Get the minimum buffer
    min_buffer = np.amin(z_buffers)
    
    total_counts = np.zeros((num_z_bins,num_z_buffers),dtype=int)
    contamination_counts = np.zeros((num_z_bins,num_z_buffers),dtype=float)
    
    # Loop over each galaxy
    for gal_i in xrange(np.shape(spec_and_phot_zs)[0]):
        
        z_s = spec_and_phot_zs[gal_i,0]
        z_p = spec_and_phot_zs[gal_i,1]
        
        # Loop over each possible lens redshift bin
        for lens_z_i in xrange(num_z_bins):
            
            lens_z = z_bins_mids[lens_z_i]
            
            # Check this galaxy's z_p is far enough away to be in ANY bin
            if z_p < lens_z + min_buffer:
                continue
            
            # Loop over buffer values
            for z_buffer_i in xrange(num_z_buffers):
                
                # Check this galaxy is far enough away
                if z_p < lens_z + z_buffers[z_buffer_i]:
                    continue
                
                # Increment the total counts for this bin and buffer
                total_counts[lens_z_i,z_buffer_i] += 1
                
                # Check if this is a contaminating galaxy
                source_z_i = bf.get_bin_index(z_s, z_bins)
                
                if(source_z_i==lens_z_i):
                    contamination_counts[lens_z_i,z_buffer_i] += 1
                    
    # Get the contamination fraction
    contamination_fractions_nd = contamination_counts/total_counts
    
    # Get the error on the contamination fraction    
    contamination_fraction_errors_nd = contamination_fractions_nd * np.sqrt( 1/contamination_counts + 1/total_counts )
    
    contamination_fractions = np.transpose(contamination_fractions_nd)
    contamination_fraction_errors = np.transpose(contamination_fraction_errors_nd)
                
    # Plot the results
    
    # Setup the fiture
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("Lens Redshift Bin Mid")
    ax.set_ylabel("Contamination Fraction")
    
    z_buffer_mean = np.mean(z_buffers)
    
    for fracs, errs, z_buffer in zip(contamination_fractions,contamination_fraction_errors,z_buffers):
        ax.errorbar(z_bins_mids+(z_buffer-z_buffer_mean)/8, fracs, yerr=errs, label=r"$\Delta z$ = "+str(z_buffer))
        
    ax.legend(numpoints=1,loc='lower right')
    ax.semilogy()
    ax.set_xlim(0.2,1.3)
    ax.set_ylim(0.0005,0.3)
    
    # Save the figure
    outfile_name = os.path.splitext(redshift_table_name)[0] + "_contamination_fraction.eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)