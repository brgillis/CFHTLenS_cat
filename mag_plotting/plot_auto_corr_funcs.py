#!/usr/bin/env python

import os
import subprocess as sbp
import sys

import astropy.io.ascii as ascii

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import bins_funcs as bf

def main(argv):
    
    # Magic values
    
    xi_min = -1.
    xi_max = 1.
    
    default_corr_funcs_table_name = "/home/brg/git/CFHTLenS_cat/Data/auto_corr_funcs.dat"
    
    default_R_col_name = "R_bin_mid_kpc"
    
    default_z_col_name = "z_bin_mid"
    default_z_bins_min = 0.2
    default_z_bins_max = 0.7
    default_num_z_bins = 5
    default_z_bins_log = False
    
    default_m_col_name = "m_bin_mid_Msun"
    default_m_bins_min = 1e9
    default_m_bins_max = 1e12
    default_num_m_bins = 3
    default_m_bins_log = True
    
    default_mag_col_name = "mag_bin_mid"
    default_mag_bins_min = 15
    default_mag_bins_max = 20
    default_num_mag_bins = 1
    default_mag_bins_log = False
    
    default_xi_col_name = "xi_mp"
    
    default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    # Load in values from command line if provided, otherwise use defaults
    
    num_args = len(argv)
    cur_arg = 2
    
    if(num_args) < cur_arg:
        corr_funcs_table_name = default_corr_funcs_table_name
    else:
        corr_funcs_table_name = argv[cur_arg]
    
    # Load in the redshift table
    try:
        corr_funcs_table = ascii.read(corr_funcs_table_name)
    except:
        print("ERROR: Corr funcs table " + corr_funcs_table + " cannot be read.")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        R_col_name = default_R_col_name
    else:
        R_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        Rs = corr_funcs_table[R_col_name]
    except:
        print("ERROR: Cannot read column " + R_col_name + " from table " + corr_funcs_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_col_name = default_z_col_name
    else:
        z_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        redshifts = corr_funcs_table[z_col_name]
    except:
        print("ERROR: Cannot read column " + z_col_name + " from table " + corr_funcs_table_name + ".")
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
        z_bins_log = default_z_bins_log
    else:
        z_bins_log = bool(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        m_col_name = default_m_col_name
    else:
        m_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        masses = corr_funcs_table[m_col_name]
    except:
        print("ERROR: Cannot read column " + m_col_name + " from table " + corr_funcs_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        m_bins_min = default_m_bins_min
    else:
        m_bins_min = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        m_bins_max = default_m_bins_max
    else:
        m_bins_max = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        num_m_bins = default_num_m_bins
    else:
        num_m_bins = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        m_bins_log = default_m_bins_log
    else:
        m_bins_log = bool(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        mag_col_name = default_mag_col_name
    else:
        mag_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        _mags = corr_funcs_table[mag_col_name]
    except:
        print("ERROR: Cannot read column " + mag_col_name + " from table " + corr_funcs_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        mag_bins_min = default_mag_bins_min
    else:
        mag_bins_min = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        mag_bins_max = default_mag_bins_max
    else:
        mag_bins_max = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        num_mag_bins = default_num_mag_bins
    else:
        num_mag_bins = float(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        mag_bins_log = default_mag_bins_log
    else:
        mag_bins_log = bool(argv[cur_arg])
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        xi_col_name = default_xi_col_name
    else:
        xi_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        xis = corr_funcs_table[xi_col_name]
    except:
        print("ERROR: Cannot read column " + xi_col_name + " from table " + corr_funcs_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        paper_location = default_paper_location
    else:
        paper_location = argv[cur_arg]
    
    # Set up bins
    if(z_bins_log):
        z_bins, z_bins_mids = bf.setup_log_bins(z_bins_min, z_bins_max, num_z_bins)
    else:
        z_bins, z_bins_mids = bf.setup_lin_bins(z_bins_min, z_bins_max, num_z_bins)
        
    if(m_bins_log):
        m_bins, m_bins_mids = bf.setup_log_bins(m_bins_min, m_bins_max, num_m_bins)
    else:
        m_bins, m_bins_mids = bf.setup_lin_bins(m_bins_min, m_bins_max, num_m_bins)
        
    if(mag_bins_log):
        _mag_bins, _mag_bins_mids = bf.setup_log_bins(mag_bins_min, mag_bins_max, num_mag_bins)
    else:
        _mag_bins, _mag_bins_mids = bf.setup_lin_bins(mag_bins_min, mag_bins_max, num_mag_bins)
    
    # Bin all values in the table
    binned_Rs = []
    binned_xis = []
    
    for z_i in xrange(num_z_bins):
        zR_list = []        
        zxi_list = []        
        for m_i in xrange(num_m_bins):
            zR_list.append([])
            zxi_list.append([])
        binned_Rs.append(zR_list)
        binned_xis.append(zxi_list)
    
    for z, m, R, xi in zip(redshifts, masses, Rs, xis):
        
        z_i = bf.get_bin_index(z, z_bins)
        if(z_i<0): continue
        
        m_i = bf.get_bin_index(m, m_bins)
        if(m_i<0): continue
        
        if(xi>xi_max): continue
        if(xi<xi_min): continue
        
        binned_Rs[z_i][m_i].append(R)
        binned_xis[z_i][m_i].append(xi)
        
    # Now start setting up the plot
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("Projected Separation (kpc)",labelpad=20)
    ax.set_ylabel(r"$\xi$",labelpad=35)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    for z_i in xrange(num_z_bins):
        
        z = z_bins_mids[z_i]
        
        for m_i in xrange(num_m_bins):
        
            m = m_bins_mids[m_i]
    
            ax = fig.add_subplot( num_m_bins, num_z_bins, z_i + num_z_bins*m_i + 1)
            ax.semilogy( binned_Rs[z_i][m_i], binned_xis[z_i][m_i], 'r' )
            #ax.legend( [emptiness,emptiness] , [r"$z_{mid}$="+ str(z) ,r"$M_{mid}$=" + "%.1E" % m],loc='upper right')
            ax.set_ylim(0.0001,0.01)
            
            # Label the redshift and mass
            xmin = 0.
            xmax = 1.
            ymin = 0.
            ymax = 1.
            ax.text(xmin+(xmax-xmin)*0.9, ymin+(ymax-ymin)*0.85, r"$z_{mid}$="+ str(z), size=10,
                    horizontalalignment='right', transform = ax.transAxes)
            ax.text(xmin+(xmax-xmin)*0.9, ymin+(ymax-ymin)*0.75, r"$M_{mid}$=" + "%.1E" % m, size=10,
                    horizontalalignment='right', transform = ax.transAxes)
            
            # set the labels as appropriate
            if(z_i!=0): # Not on the left column
                ax.set_yticklabels([])
                
            if(m_i!=num_m_bins-1): # Not on the bottom row
                ax.set_xticklabels([])
                
            if((z_i==0) and (m_i==num_m_bins-1)): # bottom-left
                ax.set_yticks([0.0001, 0.001, 0.01])
                ax.set_yticklabels([0.0001, 0.001, 0.01],fontsize=10)
                ax.set_xticks([0, 500,1000,1500,2000])
                ax.set_xticklabels([0, 500,1000,1500,2000],fontsize=10)
                continue
                
            if(z_i==0): # left column
                ax.set_yticks([0.001, 0.01])
                ax.set_yticklabels([0.001, 0.01],fontsize=10)
                
            if(m_i==num_m_bins-1): # bottom row
                ax.set_xticks([500,1000,1500,2000])
                ax.set_xticklabels([500,1000,1500,2000],fontsize=10)
    
    
    # Save the figure
    outfile_name = os.path.splitext(corr_funcs_table_name)[0] + ".eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)