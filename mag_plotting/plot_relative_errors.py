#!/usr/bin/env python

import os
import subprocess as sbp
import sys

import astropy.io.ascii as ascii

import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import bins_funcs as bf

m_sun = 1.9891e30 # kg
pc = 149597870700. * 648000. / np.pi # m

def main(argv):
    
    # Magic values
    
    figsize = (12,5)
    labelsize = 10
    
    z_shift = 0.01
    
    fitting_results_table_name = "/home/brg/git/CFHTLenS_cat/Data/gg_lensing_signal_20_bins_fitting_results.dat"
    
    paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    z_col_name = "z_min"
    z_mid_offset = 0.05
    z_mid_factor = 1
    
    m_col_name = "m_min"
    m_mid_offset = 0
    m_mid_factor = 5.5/m_sun
    
    m_bins_min = 1e9
    m_bins_max = 1e12
    num_m_bins = 3
    m_bins_log = True
    
    shear_sat_m_err_col_name = "shear_sat_m_err"
    shear_sat_m_err_offset = 0
    shear_sat_m_err_factor = 1
    
    shear_group_m_err_col_name = "shear_group_m_err"
    shear_group_m_err_offset = 0
    shear_group_m_err_factor = 1
    
    shear_sat_frac_err_col_name = "shear_sat_frac_err"
    shear_sat_frac_err_offset = 0
    shear_sat_frac_err_factor = 1
    
    magf_sat_m_err_col_name = "magf_sat_m_err"
    magf_sat_m_err_offset = 0
    magf_sat_m_err_factor = 1
    
    magf_group_m_err_col_name = "magf_group_m_err"
    magf_group_m_err_offset = 0
    magf_group_m_err_factor = 1
    
    magf_sat_frac_err_col_name = "magf_sat_frac_err"
    magf_sat_frac_err_offset = 0
    magf_sat_frac_err_factor = 1
    
    magf_Sigma_offset_err_col_name = "magf_Sigma_offset_err"
    magf_Sigma_offset_err_offset = 0
    magf_Sigma_offset_err_factor = pc*pc/m_sun
    
    overall_sat_m_err_col_name = "overall_sat_m_err"
    overall_sat_m_err_offset = 0
    overall_sat_m_err_factor = 1
    
    overall_group_m_err_col_name = "overall_group_m_err"
    overall_group_m_err_offset = 0
    overall_group_m_err_factor = 1
    
    overall_sat_frac_err_col_name = "overall_sat_frac_err"
    overall_sat_frac_err_offset = 0
    overall_sat_frac_err_factor = 1
    
    overall_Sigma_offset_err_col_name = "overall_Sigma_offset_err"
    overall_Sigma_offset_err_offset = 0
    overall_Sigma_offset_err_factor = pc*pc/m_sun
    
    
    
    # Load in the fitting results table
    try:
        fitting_results_table = ascii.read(fitting_results_table_name)
    except:
        print("ERROR: Table " + fitting_results_table_name + " cannot be read.")
        return
    
    zs = fitting_results_table[z_col_name]
    zs = np.array(zs)
    zs = (zs + z_mid_offset)*z_mid_factor
    
    ms = fitting_results_table[m_col_name]
    ms = np.array(ms)
    ms = (ms + m_mid_offset)*m_mid_factor
    
    shear_sat_m_errs = fitting_results_table[shear_sat_m_err_col_name]
    shear_sat_m_errs = np.array(shear_sat_m_errs)
    shear_sat_m_errs = (shear_sat_m_errs + shear_sat_m_err_offset)*shear_sat_m_err_factor
    
    shear_group_m_errs = fitting_results_table[shear_group_m_err_col_name]
    shear_group_m_errs = np.array(shear_group_m_errs)
    shear_group_m_errs = (shear_group_m_errs + shear_group_m_err_offset)*shear_group_m_err_factor
    
    shear_sat_frac_errs = fitting_results_table[shear_sat_frac_err_col_name]
    shear_sat_frac_errs = np.array(shear_sat_frac_errs)
    shear_sat_frac_errs = (shear_sat_frac_errs + shear_sat_frac_err_offset)*shear_sat_frac_err_factor
    
    magf_sat_m_errs = fitting_results_table[magf_sat_m_err_col_name]
    magf_sat_m_errs = np.array(magf_sat_m_errs)
    magf_sat_m_errs = (magf_sat_m_errs + magf_sat_m_err_offset)*magf_sat_m_err_factor
    
    magf_group_m_errs = fitting_results_table[magf_group_m_err_col_name]
    magf_group_m_errs = np.array(magf_group_m_errs)
    magf_group_m_errs = (magf_group_m_errs + magf_group_m_err_offset)*magf_group_m_err_factor
    
    magf_sat_frac_errs = fitting_results_table[magf_sat_frac_err_col_name]
    magf_sat_frac_errs = np.array(magf_sat_frac_errs)
    magf_sat_frac_errs = (magf_sat_frac_errs + magf_sat_frac_err_offset)*magf_sat_frac_err_factor
    
    magf_Sigma_offset_errs = fitting_results_table[magf_Sigma_offset_err_col_name]
    magf_Sigma_offset_errs = np.array(magf_Sigma_offset_errs)
    magf_Sigma_offset_errs = (magf_Sigma_offset_errs + magf_Sigma_offset_err_offset)*magf_Sigma_offset_err_factor
    
    overall_sat_m_errs = fitting_results_table[overall_sat_m_err_col_name]
    overall_sat_m_errs = np.array(overall_sat_m_errs)
    overall_sat_m_errs = (overall_sat_m_errs + overall_sat_m_err_offset)*overall_sat_m_err_factor
    
    overall_group_m_errs = fitting_results_table[overall_group_m_err_col_name]
    overall_group_m_errs = np.array(overall_group_m_errs)
    overall_group_m_errs = (overall_group_m_errs + overall_group_m_err_offset)*overall_group_m_err_factor
    
    overall_sat_frac_errs = fitting_results_table[overall_sat_frac_err_col_name]
    overall_sat_frac_errs = np.array(overall_sat_frac_errs)
    overall_sat_frac_errs = (overall_sat_frac_errs + overall_sat_frac_err_offset)*overall_sat_frac_err_factor
    
    overall_Sigma_offset_errs = fitting_results_table[overall_Sigma_offset_err_col_name]
    overall_Sigma_offset_errs = np.array(overall_Sigma_offset_errs)
    overall_Sigma_offset_errs = (overall_Sigma_offset_errs + overall_Sigma_offset_err_offset)*overall_Sigma_offset_err_factor
    
    # Set up bins
    if(m_bins_log):
        m_bins, m_bins_mids = bf.setup_log_bins(m_bins_min, m_bins_max, num_m_bins)
    else:
        m_bins, m_bins_mids = bf.setup_lin_bins(m_bins_min, m_bins_max, num_m_bins)
    
    # Bin all values in the table
    binned_zs = []
    binned_shear_sat_m_errs = []
    binned_shear_group_m_errs = []
    binned_shear_sat_frac_errs = []
    binned_magf_sat_m_errs = []
    binned_magf_group_m_errs = []
    binned_magf_sat_frac_errs = []
    binned_magf_Sigma_offset_errs = []
    binned_overall_sat_m_errs = []
    binned_overall_group_m_errs = []
    binned_overall_sat_frac_errs = []
    binned_overall_Sigma_offset_errs = []
    
    for m_i in xrange(num_m_bins):
        binned_zs.append([])
        binned_shear_sat_m_errs.append([])
        binned_shear_group_m_errs.append([])
        binned_shear_sat_frac_errs.append([])
        binned_magf_sat_m_errs.append([])
        binned_magf_group_m_errs.append([])
        binned_magf_sat_frac_errs.append([])
        binned_magf_Sigma_offset_errs.append([])
        binned_overall_sat_m_errs.append([])
        binned_overall_group_m_errs.append([])
        binned_overall_sat_frac_errs.append([])
        binned_overall_Sigma_offset_errs.append([])

    
    for z, m, shear_sat_m_err, shear_group_m_err, shear_sat_frac_err, \
        magf_sat_m_err, magf_group_m_err, magf_sat_frac_err, magf_Sigma_offset_err, \
        overall_sat_m_err, overall_group_m_err, overall_sat_frac_err, overall_Sigma_offset_err \
         in zip(zs, ms, shear_sat_m_errs, shear_group_m_errs, shear_sat_frac_errs, \
                magf_sat_m_errs, magf_group_m_errs, magf_sat_frac_errs, magf_Sigma_offset_errs, \
                overall_sat_m_errs, overall_group_m_errs, overall_sat_frac_errs, overall_Sigma_offset_errs):
        
        m_i = bf.get_bin_index(m, m_bins)
        if(m_i<0): continue
        
        binned_zs[m_i].append(z)
        binned_shear_sat_m_errs[m_i].append(shear_sat_m_err)
        binned_shear_group_m_errs[m_i].append(shear_group_m_err)
        binned_shear_sat_frac_errs[m_i].append(shear_sat_frac_err)
        binned_magf_sat_m_errs[m_i].append(magf_sat_m_err)
        binned_magf_group_m_errs[m_i].append(magf_group_m_err)
        binned_magf_sat_frac_errs[m_i].append(magf_sat_frac_err)
        binned_magf_Sigma_offset_errs[m_i].append(magf_Sigma_offset_err)
        binned_overall_sat_m_errs[m_i].append(overall_sat_m_err)
        binned_overall_group_m_errs[m_i].append(overall_group_m_err)
        binned_overall_sat_frac_errs[m_i].append(overall_sat_frac_err)
        binned_overall_Sigma_offset_errs[m_i].append(overall_Sigma_offset_err)
        
    # Convert all lists to arrays      
    for m_i in xrange(num_m_bins):
        binned_zs= np.array(binned_zs)
        binned_shear_sat_m_errs= np.array(binned_shear_sat_m_errs)
        binned_shear_group_m_errs= np.array(binned_shear_group_m_errs)
        binned_shear_sat_frac_errs= np.array(binned_shear_sat_frac_errs)
        binned_magf_sat_m_errs= np.array(binned_magf_sat_m_errs)
        binned_magf_group_m_errs= np.array(binned_magf_group_m_errs)
        binned_magf_sat_frac_errs= np.array(binned_magf_sat_frac_errs)
        binned_magf_Sigma_offset_errs= np.array(binned_magf_Sigma_offset_errs)
        binned_overall_sat_m_errs= np.array(binned_overall_sat_m_errs)
        binned_overall_group_m_errs= np.array(binned_overall_group_m_errs)
        binned_overall_sat_frac_errs= np.array(binned_overall_sat_frac_errs)
        binned_overall_Sigma_offset_errs= np.array(binned_overall_Sigma_offset_errs)
        
    # Do the shear plot now
    fig = pyplot.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")
    ax.set_xlabel("Redshift",labelpad=20)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    num_p_bins = 4
        
    for m_i in xrange(num_m_bins):
    
        m = m_bins_mids[m_i]
        
        for p_i in xrange(num_p_bins):
            
            ax = fig.add_subplot( num_m_bins, num_p_bins, p_i + num_p_bins*m_i + 1)
        
            ax.set_xlim(0.2,1.3)
            
            if(p_i==0):
                # Plot sat_m
                
                ax.plot([0.2,1.3],[1,1],label=None,color="b",linestyle="solid")
                
                ax.plot( binned_zs[m_i], binned_magf_sat_m_errs[m_i]/binned_shear_sat_m_errs[m_i],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r')
                
                ideal_err = ( binned_magf_sat_m_errs[m_i]**-2 + binned_shear_sat_m_errs[m_i]**-2 )**-0.5
                
                ax.plot( binned_zs[m_i], ideal_err/binned_shear_sat_m_errs[m_i],
                             color='m', linestyle='None', label="Ideal errors",
                             marker='x')
                
                ax.plot( binned_zs[m_i]+z_shift, binned_overall_sat_m_errs[m_i]/binned_shear_sat_m_errs[m_i],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*')
                
                ax.set_ylabel(r"$\sigma(M_{\rm 1h})/\sigma_{\rm shear}(M_{\rm 1h})$",labelpad=5)
            elif(p_i==1):
                # Plot group_m
                ax.plot([0.2,1.3],[1,1],label=None,color="b",linestyle="solid")
                
                ax.plot( binned_zs[m_i], binned_magf_group_m_errs[m_i]/binned_shear_group_m_errs[m_i],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r')
                
                ideal_err = ( binned_magf_group_m_errs[m_i]**-2 + binned_shear_group_m_errs[m_i]**-2 )**-0.5
                
                ax.plot( binned_zs[m_i], ideal_err/binned_shear_group_m_errs[m_i],
                             color='m', linestyle='None', label="Ideal errors",
                             marker='x')
                
                ax.plot( binned_zs[m_i]+z_shift, binned_overall_group_m_errs[m_i]/binned_shear_group_m_errs[m_i],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*')
                
                ax.set_ylabel(r"$\sigma(M_{\rm gr})/\sigma_{\rm shear}(M_{\rm gr})$",labelpad=5)
            elif(p_i==2):
                # Plot sat_frac
                ax.plot([0.2,1.3],[1,1],label=None,color="b",linestyle="solid")
                ax.plot( binned_zs[m_i], binned_magf_sat_frac_errs[m_i]/binned_shear_sat_frac_errs[m_i],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r')
                
                ideal_err = ( binned_magf_sat_frac_errs[m_i]**-2 + binned_shear_sat_frac_errs[m_i]**-2 )**-0.5
                
                ax.plot( binned_zs[m_i], ideal_err/binned_shear_sat_frac_errs[m_i],
                             color='m', linestyle='None', label="Ideal errors",
                             marker='x')
                
                
                ax.plot( binned_zs[m_i]+z_shift, binned_overall_sat_frac_errs[m_i]/binned_shear_sat_frac_errs[m_i],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*')
                ax.set_ylabel(r"$\sigma(f_{\rm sat})/\sigma_{\rm shear}(f_{\rm sat})$",labelpad=5)
            elif(p_i==3):
                # Plot product
                shear_prod = binned_shear_sat_m_errs[m_i]*binned_shear_group_m_errs[m_i]*binned_shear_sat_frac_errs[m_i]
                magf_prod = binned_magf_sat_m_errs[m_i]*binned_magf_group_m_errs[m_i]*binned_magf_sat_frac_errs[m_i]
                
                ax.plot([0.2,1.3],[1,1],label=None,color="b",linestyle="solid")
                
                ax.plot( binned_zs[m_i],  magf_prod/shear_prod,
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r')
                
                ideal_prod = ( magf_prod**-2 + shear_prod**-2 )**-0.5
                ax.plot( binned_zs[m_i],  ideal_prod/shear_prod,
                             color='m', linestyle='None', label="Ideal errors",
                             marker='x')
                
                ax.plot( binned_zs[m_i],
                         binned_overall_sat_m_errs[m_i]*binned_overall_group_m_errs[m_i]*binned_overall_sat_frac_errs[m_i]/shear_prod,
                             color='k', linestyle='None', label="Overall fit",
                             marker='*')
                ax.set_ylabel(r"$\sigma^3/\sigma_{\rm shear}^3$",labelpad=5)
            
            if(p_i==3):
                # Special y-axis
                ax.set_ylim(0.05**3,10**3)
                ax.set_yscale("log", nonposy='clip')
                ax.set_yticks(np.array([0.1, 0.3, 1.0, 3.0, 10])**3)
                ax.set_yticklabels(np.array([0.1, 0.3, 1.0, 3.0, 10])**3,fontsize=12)
                
            else:
                ax.set_ylim(0.05,10)
                ax.set_yscale("log", nonposy='clip')
                ax.set_yticks([0.1, 0.3, 1.0, 3.0, 10])
                ax.set_yticklabels([0.1, 0.3, 1.0, 3.0, 10],fontsize=12)
            
            # Label the mass
            xmin = 0.
            xmax = 1.
            ymin = 0.
            ymax = 1.
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.9, r"$M_{mid}$=" + "%.1e" % m, size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            
            # set the x ticks and labels as appropriate
            ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
                            
            if(m_i!=num_m_bins-1): # Not on the bottom row
                ax.set_xticklabels([])
                
            if(m_i==num_m_bins-1): # bottom row
                ax.set_xticklabels([.2, .4, .6, .8, 1.0, 1.2],fontsize=12)
    
    
    # Save the figure
    outfile_name = os.path.splitext(fitting_results_table_name)[0] + "_relative_errors.eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)