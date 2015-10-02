#!/usr/bin/env python

import matplotlib
import os
import sys

import icebrgpy.bins_funcs as bf
import icebrgpy.columns_reader as cr
import matplotlib.pyplot as pyplot
import numpy as np
import subprocess as sbp


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


m_sun = 1.9891e30 # kg
pc = 149597870700. * 648000. / np.pi # m

sqdeg_per_sqrad = (180/np.pi)**2

def main(argv):
    
    # Magic values
    
    unweighted_table_name = "/home/brg/git/CFHTLenS_cat/Data/gg_lensing_signal_20_bins_unweighted_fitting_results.dat"
    weighted_table_name = "/home/brg/git/CFHTLenS_cat/Data/gg_lensing_signal_20_bins_fitting_results.dat"
    
    paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    m_bins_min = 1e9
    m_bins_max = 1e12
    num_m_bins = 3
    m_bins_log = True
    
    unweighted_reader = cr.columns_reader()
    weighted_reader = cr.columns_reader()
    
    for reader in [unweighted_reader, weighted_reader]:
    
        reader.add("z","z_min",offset=0.05)
        reader.add("m","m_min",factor=1./m_sun)
        reader.add("sigma_offset_best","magf_Sigma_offset_best",factor=pc*pc/m_sun)
        reader.add("sigma_offset_err","magf_Sigma_offset_err",factor=pc*pc/m_sun)
        reader.add("sigma_crit","Sigma_crit",factor=pc*pc/m_sun)
        
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("Lens Redshift Bin Centre",fontsize=16,labelpad=15)
    
    ax.set_ylabel(r"Change in $\kappa_{\rm offset}$ from field weighting",fontsize=16,labelpad=40)
    
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    unweighted_cols = unweighted_reader.read(unweighted_table_name)
    weighted_cols = weighted_reader.read(weighted_table_name)
    
    binned_uw_cols = []
    binned_w_cols = []
        
    for col_i in xrange(unweighted_reader.num_cols()):
        binned_uw_cols.append([])
        binned_w_cols.append([])
        
    
    # Set up bins
    if(m_bins_log):
        m_bins, m_bins_mids = bf.setup_log_bins(m_bins_min, m_bins_max, num_m_bins)
    else:
        m_bins, m_bins_mids = bf.setup_lin_bins(m_bins_min, m_bins_max, num_m_bins)
    
    # Bin all values in the table
    for col_i in xrange(unweighted_reader.num_cols()):
        for _ in xrange(num_m_bins):
            binned_uw_cols[col_i].append([])
            binned_w_cols[col_i].append([])

    for elems in zip(*unweighted_cols):
        
        m = elems[unweighted_reader.index("m")]
        
        m_i = bf.get_bin_index(m+0.001, m_bins)
        if(m_i<0): continue
        
        for col_i in xrange(unweighted_reader.num_cols()):
            binned_uw_cols[col_i][m_i].append(elems[col_i])

    for elems in zip(*weighted_cols):
        
        m = elems[weighted_reader.index("m")]
        
        m_i = bf.get_bin_index(m+0.001, m_bins)
        if(m_i<0): continue
        
        for col_i in xrange(unweighted_reader.num_cols()):
            binned_w_cols[col_i][m_i].append(elems[col_i])
        
    # Convert all lists to arrays      
    for col_i in xrange(unweighted_reader.num_cols()):
        for m_i in xrange(num_m_bins):
            binned_uw_cols[col_i][m_i] = np.array(binned_uw_cols[col_i][m_i])
            binned_w_cols[col_i][m_i] = np.array(binned_w_cols[col_i][m_i])
        
    # Do the plotting now
    
    for m_i in xrange(num_m_bins):
        m = m_bins_mids[m_i]
        
        ax = fig.add_subplot( num_m_bins, 1, m_i+1)
        
        sigma_off_diff = (binned_w_cols[weighted_reader.index("sigma_offset_best")][m_i] -
                          binned_uw_cols[unweighted_reader.index("sigma_offset_best")][m_i])
        sigma_off_err = np.sqrt(binned_w_cols[weighted_reader.index("sigma_offset_err")][m_i]**2 +
                         binned_uw_cols[unweighted_reader.index("sigma_offset_err")][m_i]**2)
        
        kappa_off_diff = sigma_off_diff/binned_w_cols[weighted_reader.index("sigma_crit")][m_i]
        kappa_off_err = sigma_off_err/binned_w_cols[weighted_reader.index("sigma_crit")][m_i]
                         
        ax.errorbar( binned_w_cols[unweighted_reader.index("z")][m_i], kappa_off_diff,
                             linestyle='None', color='r',
                             marker='o', markerfacecolor='none', markeredgecolor='r',yerr=kappa_off_err )
        ax.plot([0.2,1.3],[0,0],label=None,color="k",linestyle="dashed")
        ax.set_ylim(-0.00549,0.00549)
        ax.set_xlim(0.2,1.3)
            
        # Label the mass
        xmin = 0.
        xmax = 1.
        ymin = 0.
        ymax = 1.
        ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.9, r"$M_{\rm mid}$=" + "%.1e" % m, size=12,
                horizontalalignment='right', transform = ax.transAxes)
                            
        if(m_i!=num_m_bins-1): # Not on the bottom row
            ax.set_xticklabels([])
        
        
    # Save the figure
    outfile_name = os.path.splitext(unweighted_table_name)[0] + "_weight_comp.eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
            
    # Show the figure
        
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)