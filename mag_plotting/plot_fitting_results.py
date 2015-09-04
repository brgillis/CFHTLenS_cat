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
    
    figsize = (12,5)
    labelsize = 10
    
    z_shift = 0.01
    
    chi_2_threshold = 30
    overall_chi_2_threshold = 50
    
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
    
    shear_chi_squared_col_name = "shear_Chi_squared"
    shear_chi_squared_offset = 0
    shear_chi_squared_factor = 1
    
    shear_sat_m_best_col_name = "shear_sat_m_best"
    shear_sat_m_best_offset = -np.log10(m_sun)
    shear_sat_m_best_factor = 1
    
    shear_sat_m_err_col_name = "shear_sat_m_err"
    shear_sat_m_err_offset = 0
    shear_sat_m_err_factor = 1
    
    shear_group_m_best_col_name = "shear_group_m_best"
    shear_group_m_best_offset = -np.log10(m_sun)
    shear_group_m_best_factor = 1
    
    shear_group_m_err_col_name = "shear_group_m_err"
    shear_group_m_err_offset = 0
    shear_group_m_err_factor = 1
    
    shear_sat_frac_best_col_name = "shear_sat_frac_best"
    shear_sat_frac_best_offset = 0
    shear_sat_frac_best_factor = 1
    
    shear_sat_frac_err_col_name = "shear_sat_frac_err"
    shear_sat_frac_err_offset = 0
    shear_sat_frac_err_factor = 1
    
    magf_chi_squared_col_name = "magf_Chi_squared"
    magf_chi_squared_offset = 0
    magf_chi_squared_factor = 1
    
    magf_sat_m_best_col_name = "magf_sat_m_best"
    magf_sat_m_best_offset = -np.log10(m_sun)
    magf_sat_m_best_factor = 1
    
    magf_sat_m_err_col_name = "magf_sat_m_err"
    magf_sat_m_err_offset = 0
    magf_sat_m_err_factor = 1
    
    magf_group_m_best_col_name = "magf_group_m_best"
    magf_group_m_best_offset = -np.log10(m_sun)
    magf_group_m_best_factor = 1
    
    magf_group_m_err_col_name = "magf_group_m_err"
    magf_group_m_err_offset = 0
    magf_group_m_err_factor = 1
    
    magf_sat_frac_best_col_name = "magf_sat_frac_best"
    magf_sat_frac_best_offset = 0
    magf_sat_frac_best_factor = 1
    
    magf_sat_frac_err_col_name = "magf_sat_frac_err"
    magf_sat_frac_err_offset = 0
    magf_sat_frac_err_factor = 1
    
    magf_Sigma_offset_best_col_name = "magf_Sigma_offset_best"
    magf_Sigma_offset_best_offset = 0
    magf_Sigma_offset_best_factor = pc*pc/m_sun
    
    magf_Sigma_offset_err_col_name = "magf_Sigma_offset_err"
    magf_Sigma_offset_err_offset = 0
    magf_Sigma_offset_err_factor = pc*pc/m_sun
    
    overall_chi_squared_col_name = "overall_Chi_squared"
    overall_chi_squared_offset = 0
    overall_chi_squared_factor = 0.5
    
    overall_sat_m_best_col_name = "overall_sat_m_best"
    overall_sat_m_best_offset = -np.log10(m_sun)
    overall_sat_m_best_factor = 1
    
    overall_sat_m_err_col_name = "overall_sat_m_err"
    overall_sat_m_err_offset = 0
    overall_sat_m_err_factor = 1
    
    overall_group_m_best_col_name = "overall_group_m_best"
    overall_group_m_best_offset = -np.log10(m_sun)
    overall_group_m_best_factor = 1
    
    overall_group_m_err_col_name = "overall_group_m_err"
    overall_group_m_err_offset = 0
    overall_group_m_err_factor = 1
    
    overall_sat_frac_best_col_name = "overall_sat_frac_best"
    overall_sat_frac_best_offset = 0
    overall_sat_frac_best_factor = 1
    
    overall_sat_frac_err_col_name = "overall_sat_frac_err"
    overall_sat_frac_err_offset = 0
    overall_sat_frac_err_factor = 1
    
    overall_Sigma_offset_best_col_name = "overall_Sigma_offset_best"
    overall_Sigma_offset_best_offset = 0
    overall_Sigma_offset_best_factor = pc*pc/m_sun
    
    overall_Sigma_offset_err_col_name = "overall_Sigma_offset_err"
    overall_Sigma_offset_err_offset = 0
    overall_Sigma_offset_err_factor = pc*pc/m_sun
    
    Sigma_crit_col_name = "Sigma_crit"
    Sigma_crit_offset = 0
    Sigma_crit_factor = pc*pc/m_sun
    
    
    
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
    
    shear_chi_squareds = fitting_results_table[shear_chi_squared_col_name]
    shear_chi_squareds = np.array(shear_chi_squareds)
    shear_chi_squareds = (shear_chi_squareds + shear_chi_squared_offset)*shear_chi_squared_factor
    
    shear_sat_m_bests = fitting_results_table[shear_sat_m_best_col_name]
    shear_sat_m_bests = np.array(shear_sat_m_bests)
    shear_sat_m_bests = (shear_sat_m_bests + shear_sat_m_best_offset)*shear_sat_m_best_factor
    
    shear_sat_m_errs = fitting_results_table[shear_sat_m_err_col_name]
    shear_sat_m_errs = np.array(shear_sat_m_errs)
    shear_sat_m_errs = (shear_sat_m_errs + shear_sat_m_err_offset)*shear_sat_m_err_factor
    
    shear_group_m_bests = fitting_results_table[shear_group_m_best_col_name]
    shear_group_m_bests = np.array(shear_group_m_bests)
    shear_group_m_bests = (shear_group_m_bests + shear_group_m_best_offset)*shear_group_m_best_factor
    
    shear_group_m_errs = fitting_results_table[shear_group_m_err_col_name]
    shear_group_m_errs = np.array(shear_group_m_errs)
    shear_group_m_errs = (shear_group_m_errs + shear_group_m_err_offset)*shear_group_m_err_factor
    
    shear_sat_frac_bests = fitting_results_table[shear_sat_frac_best_col_name]
    shear_sat_frac_bests = np.array(shear_sat_frac_bests)
    shear_sat_frac_bests = (shear_sat_frac_bests + shear_sat_frac_best_offset)*shear_sat_frac_best_factor
    
    shear_sat_frac_errs = fitting_results_table[shear_sat_frac_err_col_name]
    shear_sat_frac_errs = np.array(shear_sat_frac_errs)
    shear_sat_frac_errs = (shear_sat_frac_errs + shear_sat_frac_err_offset)*shear_sat_frac_err_factor
    
    magf_chi_squareds = fitting_results_table[magf_chi_squared_col_name]
    magf_chi_squareds = np.array(magf_chi_squareds)
    magf_chi_squareds = (magf_chi_squareds + magf_chi_squared_offset)*magf_chi_squared_factor
    
    magf_sat_m_bests = fitting_results_table[magf_sat_m_best_col_name]
    magf_sat_m_bests = np.array(magf_sat_m_bests)
    magf_sat_m_bests = (magf_sat_m_bests + magf_sat_m_best_offset)*magf_sat_m_best_factor
    
    magf_sat_m_errs = fitting_results_table[magf_sat_m_err_col_name]
    magf_sat_m_errs = np.array(magf_sat_m_errs)
    magf_sat_m_errs = (magf_sat_m_errs + magf_sat_m_err_offset)*magf_sat_m_err_factor
    
    magf_group_m_bests = fitting_results_table[magf_group_m_best_col_name]
    magf_group_m_bests = np.array(magf_group_m_bests)
    magf_group_m_bests = (magf_group_m_bests + magf_group_m_best_offset)*magf_group_m_best_factor
    
    magf_group_m_errs = fitting_results_table[magf_group_m_err_col_name]
    magf_group_m_errs = np.array(magf_group_m_errs)
    magf_group_m_errs = (magf_group_m_errs + magf_group_m_err_offset)*magf_group_m_err_factor
    
    magf_sat_frac_bests = fitting_results_table[magf_sat_frac_best_col_name]
    magf_sat_frac_bests = np.array(magf_sat_frac_bests)
    magf_sat_frac_bests = (magf_sat_frac_bests + magf_sat_frac_best_offset)*magf_sat_frac_best_factor
    
    magf_sat_frac_errs = fitting_results_table[magf_sat_frac_err_col_name]
    magf_sat_frac_errs = np.array(magf_sat_frac_errs)
    magf_sat_frac_errs = (magf_sat_frac_errs + magf_sat_frac_err_offset)*magf_sat_frac_err_factor
    
    magf_Sigma_offset_bests = fitting_results_table[magf_Sigma_offset_best_col_name]
    magf_Sigma_offset_bests = np.array(magf_Sigma_offset_bests)
    magf_Sigma_offset_bests = (magf_Sigma_offset_bests + magf_Sigma_offset_best_offset)*magf_Sigma_offset_best_factor
    
    magf_Sigma_offset_errs = fitting_results_table[magf_Sigma_offset_err_col_name]
    magf_Sigma_offset_errs = np.array(magf_Sigma_offset_errs)
    magf_Sigma_offset_errs = (magf_Sigma_offset_errs + magf_Sigma_offset_err_offset)*magf_Sigma_offset_err_factor
    
    overall_chi_squareds = fitting_results_table[overall_chi_squared_col_name]
    overall_chi_squareds = np.array(overall_chi_squareds)
    overall_chi_squareds = (overall_chi_squareds + overall_chi_squared_offset)*overall_chi_squared_factor
    
    overall_sat_m_bests = fitting_results_table[overall_sat_m_best_col_name]
    overall_sat_m_bests = np.array(overall_sat_m_bests)
    overall_sat_m_bests = (overall_sat_m_bests + overall_sat_m_best_offset)*overall_sat_m_best_factor
    
    overall_sat_m_errs = fitting_results_table[overall_sat_m_err_col_name]
    overall_sat_m_errs = np.array(overall_sat_m_errs)
    overall_sat_m_errs = (overall_sat_m_errs + overall_sat_m_err_offset)*overall_sat_m_err_factor
    
    overall_group_m_bests = fitting_results_table[overall_group_m_best_col_name]
    overall_group_m_bests = np.array(overall_group_m_bests)
    overall_group_m_bests = (overall_group_m_bests + overall_group_m_best_offset)*overall_group_m_best_factor
    
    overall_group_m_errs = fitting_results_table[overall_group_m_err_col_name]
    overall_group_m_errs = np.array(overall_group_m_errs)
    overall_group_m_errs = (overall_group_m_errs + overall_group_m_err_offset)*overall_group_m_err_factor
    
    overall_sat_frac_bests = fitting_results_table[overall_sat_frac_best_col_name]
    overall_sat_frac_bests = np.array(overall_sat_frac_bests)
    overall_sat_frac_bests = (overall_sat_frac_bests + overall_sat_frac_best_offset)*overall_sat_frac_best_factor
    
    overall_sat_frac_errs = fitting_results_table[overall_sat_frac_err_col_name]
    overall_sat_frac_errs = np.array(overall_sat_frac_errs)
    overall_sat_frac_errs = (overall_sat_frac_errs + overall_sat_frac_err_offset)*overall_sat_frac_err_factor
    
    overall_Sigma_offset_bests = fitting_results_table[overall_Sigma_offset_best_col_name]
    overall_Sigma_offset_bests = np.array(overall_Sigma_offset_bests)
    overall_Sigma_offset_bests = (overall_Sigma_offset_bests + overall_Sigma_offset_best_offset)*overall_Sigma_offset_best_factor
    
    overall_Sigma_offset_errs = fitting_results_table[overall_Sigma_offset_err_col_name]
    overall_Sigma_offset_errs = np.array(overall_Sigma_offset_errs)
    overall_Sigma_offset_errs = (overall_Sigma_offset_errs + overall_Sigma_offset_err_offset)*overall_Sigma_offset_err_factor
    
    Sigma_crits = fitting_results_table[Sigma_crit_col_name]
    Sigma_crits = np.array(Sigma_crits)
    Sigma_crits = (Sigma_crits + Sigma_crit_offset)*Sigma_crit_factor
    
    # Set up bins
    if(m_bins_log):
        m_bins, m_bins_mids = bf.setup_log_bins(m_bins_min, m_bins_max, num_m_bins)
    else:
        m_bins, m_bins_mids = bf.setup_lin_bins(m_bins_min, m_bins_max, num_m_bins)
    
    # Bin all values in the table
    binned_zs = []
    binned_shear_chi_squareds = []
    binned_shear_sat_m_bests = []
    binned_shear_sat_m_errs = []
    binned_shear_group_m_bests = []
    binned_shear_group_m_errs = []
    binned_shear_sat_frac_bests = []
    binned_shear_sat_frac_errs = []
    binned_magf_chi_squareds = []
    binned_magf_sat_m_bests = []
    binned_magf_sat_m_errs = []
    binned_magf_group_m_bests = []
    binned_magf_group_m_errs = []
    binned_magf_sat_frac_bests = []
    binned_magf_sat_frac_errs = []
    binned_magf_Sigma_offset_bests = []
    binned_magf_Sigma_offset_errs = []
    binned_overall_chi_squareds = []
    binned_overall_sat_m_bests = []
    binned_overall_sat_m_errs = []
    binned_overall_group_m_bests = []
    binned_overall_group_m_errs = []
    binned_overall_sat_frac_bests = []
    binned_overall_sat_frac_errs = []
    binned_overall_Sigma_offset_bests = []
    binned_overall_Sigma_offset_errs = []
    binned_Sigma_crits = []
    
    for m_i in xrange(num_m_bins):
        binned_zs.append([])
        binned_shear_chi_squareds.append([])
        binned_shear_sat_m_bests.append([])
        binned_shear_sat_m_errs.append([])
        binned_shear_group_m_bests.append([])
        binned_shear_group_m_errs.append([])
        binned_shear_sat_frac_bests.append([])
        binned_shear_sat_frac_errs.append([])
        binned_magf_chi_squareds.append([])
        binned_magf_sat_m_bests.append([])
        binned_magf_sat_m_errs.append([])
        binned_magf_group_m_bests.append([])
        binned_magf_group_m_errs.append([])
        binned_magf_sat_frac_bests.append([])
        binned_magf_sat_frac_errs.append([])
        binned_magf_Sigma_offset_bests.append([])
        binned_magf_Sigma_offset_errs.append([])
        binned_overall_chi_squareds.append([])
        binned_overall_sat_m_bests.append([])
        binned_overall_sat_m_errs.append([])
        binned_overall_group_m_bests.append([])
        binned_overall_group_m_errs.append([])
        binned_overall_sat_frac_bests.append([])
        binned_overall_sat_frac_errs.append([])
        binned_overall_Sigma_offset_bests.append([])
        binned_overall_Sigma_offset_errs.append([])
        binned_Sigma_crits.append([])

    
    for z, m, shear_chi_squared, shear_sat_m_best, shear_sat_m_err, shear_group_m_best, shear_group_m_err, shear_sat_frac_best, shear_sat_frac_err, \
        magf_chi_squared, magf_sat_m_best, magf_sat_m_err, magf_group_m_best, magf_group_m_err, magf_sat_frac_best, magf_sat_frac_err, magf_Sigma_offset_best, magf_Sigma_offset_err, \
        overall_chi_squared, overall_sat_m_best, overall_sat_m_err, overall_group_m_best, overall_group_m_err, overall_sat_frac_best, overall_sat_frac_err, overall_Sigma_offset_best, overall_Sigma_offset_err, \
        Sigma_crit \
         in zip(zs, ms, shear_chi_squareds, shear_sat_m_bests, shear_sat_m_errs, shear_group_m_bests, shear_group_m_errs, shear_sat_frac_bests, shear_sat_frac_errs, \
                magf_chi_squareds, magf_sat_m_bests, magf_sat_m_errs, magf_group_m_bests, magf_group_m_errs, magf_sat_frac_bests, magf_sat_frac_errs, magf_Sigma_offset_bests, magf_Sigma_offset_errs, \
                overall_chi_squareds, overall_sat_m_bests, overall_sat_m_errs, overall_group_m_bests, overall_group_m_errs, overall_sat_frac_bests, overall_sat_frac_errs, overall_Sigma_offset_bests, overall_Sigma_offset_errs, \
                Sigma_crits):
        
        m_i = bf.get_bin_index(m, m_bins)
        if(m_i<0): continue
        
        binned_zs[m_i].append(z)
        binned_shear_chi_squareds[m_i].append(shear_chi_squared)
        binned_shear_sat_m_bests[m_i].append(shear_sat_m_best)
        binned_shear_sat_m_errs[m_i].append(shear_sat_m_err)
        binned_shear_group_m_bests[m_i].append(shear_group_m_best)
        binned_shear_group_m_errs[m_i].append(shear_group_m_err)
        binned_shear_sat_frac_bests[m_i].append(shear_sat_frac_best)
        binned_shear_sat_frac_errs[m_i].append(shear_sat_frac_err)
        binned_magf_chi_squareds[m_i].append(magf_chi_squared)
        binned_magf_sat_m_bests[m_i].append(magf_sat_m_best)
        binned_magf_sat_m_errs[m_i].append(magf_sat_m_err)
        binned_magf_group_m_bests[m_i].append(magf_group_m_best)
        binned_magf_group_m_errs[m_i].append(magf_group_m_err)
        binned_magf_sat_frac_bests[m_i].append(magf_sat_frac_best)
        binned_magf_sat_frac_errs[m_i].append(magf_sat_frac_err)
        binned_magf_Sigma_offset_bests[m_i].append(magf_Sigma_offset_best)
        binned_magf_Sigma_offset_errs[m_i].append(magf_Sigma_offset_err)
        binned_overall_chi_squareds[m_i].append(overall_chi_squared)
        binned_overall_sat_m_bests[m_i].append(overall_sat_m_best)
        binned_overall_sat_m_errs[m_i].append(overall_sat_m_err)
        binned_overall_group_m_bests[m_i].append(overall_group_m_best)
        binned_overall_group_m_errs[m_i].append(overall_group_m_err)
        binned_overall_sat_frac_bests[m_i].append(overall_sat_frac_best)
        binned_overall_sat_frac_errs[m_i].append(overall_sat_frac_err)
        binned_overall_Sigma_offset_bests[m_i].append(overall_Sigma_offset_best)
        binned_overall_Sigma_offset_errs[m_i].append(overall_Sigma_offset_err)
        binned_Sigma_crits[m_i].append(Sigma_crit)
        
    # Convert all lists to arrays      
    for m_i in xrange(num_m_bins):
        binned_zs= np.array(binned_zs)
        binned_shear_chi_squareds= np.array(binned_shear_chi_squareds)
        binned_shear_sat_m_bests= np.array(binned_shear_sat_m_bests)
        binned_shear_sat_m_errs= np.array(binned_shear_sat_m_errs)
        binned_shear_group_m_bests= np.array(binned_shear_group_m_bests)
        binned_shear_group_m_errs= np.array(binned_shear_group_m_errs)
        binned_shear_sat_frac_bests= np.array(binned_shear_sat_frac_bests)
        binned_shear_sat_frac_errs= np.array(binned_shear_sat_frac_errs)
        binned_magf_chi_squareds= np.array(binned_magf_chi_squareds)
        binned_magf_sat_m_bests= np.array(binned_magf_sat_m_bests)
        binned_magf_sat_m_errs= np.array(binned_magf_sat_m_errs)
        binned_magf_group_m_bests= np.array(binned_magf_group_m_bests)
        binned_magf_group_m_errs= np.array(binned_magf_group_m_errs)
        binned_magf_sat_frac_bests= np.array(binned_magf_sat_frac_bests)
        binned_magf_sat_frac_errs= np.array(binned_magf_sat_frac_errs)
        binned_magf_Sigma_offset_bests= np.array(binned_magf_Sigma_offset_bests)
        binned_magf_Sigma_offset_errs= np.array(binned_magf_Sigma_offset_errs)
        binned_overall_chi_squareds= np.array(binned_overall_chi_squareds)
        binned_overall_sat_m_bests= np.array(binned_overall_sat_m_bests)
        binned_overall_sat_m_errs= np.array(binned_overall_sat_m_errs)
        binned_overall_group_m_bests= np.array(binned_overall_group_m_bests)
        binned_overall_group_m_errs= np.array(binned_overall_group_m_errs)
        binned_overall_sat_frac_bests= np.array(binned_overall_sat_frac_bests)
        binned_overall_sat_frac_errs= np.array(binned_overall_sat_frac_errs)
        binned_overall_Sigma_offset_bests= np.array(binned_overall_Sigma_offset_bests)
        binned_overall_Sigma_offset_errs= np.array(binned_overall_Sigma_offset_errs)
        binned_Sigma_crits= np.array(binned_Sigma_crits)
        
    # Do the plot now
    fig = pyplot.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")
    ax.set_xlabel("Lens Redshift Bin Centre",labelpad=20)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    num_p_bins = 4
        
    for m_i in xrange(num_m_bins):
    
        m = m_bins_mids[m_i]
        
        shear_good = binned_shear_chi_squareds[m_i] < chi_2_threshold
        magf_good = binned_magf_chi_squareds[m_i] < chi_2_threshold
        overall_good = binned_overall_chi_squareds[m_i] < overall_chi_2_threshold
        
        for p_i in xrange(num_p_bins):
        
            if(p_i==0):
                # Plot sat_m
                ax = fig.add_subplot( num_m_bins, num_p_bins, p_i + num_p_bins*m_i + 1)
                ax.errorbar( binned_zs[m_i][shear_good]-z_shift, binned_shear_sat_m_bests[m_i][shear_good],
                             color='b', linestyle='None', label="Shear fit",
                             marker='.', yerr=binned_shear_sat_m_errs[m_i][shear_good] )
                ax.errorbar( binned_zs[m_i][magf_good], binned_magf_sat_m_bests[m_i][magf_good],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r', markersize=3,
                             yerr=binned_magf_sat_m_errs[m_i][magf_good] )
                ax.errorbar( binned_zs[m_i][overall_good]+z_shift, binned_overall_sat_m_bests[m_i][overall_good],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*',yerr=binned_overall_sat_m_errs[m_i][overall_good] )
                ax.set_ylim(9.1,13.9)
                ax.set_ylabel(r"$\log_{10}(M_{\rm 1h}/M_{\rm sun})$",labelpad=5)
            elif(p_i==1):
                # Plot group_m
                ax = fig.add_subplot( num_m_bins, num_p_bins, p_i + num_p_bins*m_i + 1)
                ax.errorbar( binned_zs[m_i][shear_good]-z_shift, binned_shear_group_m_bests[m_i][shear_good],
                             color='b', linestyle='None', label="Shear fit",
                             marker='.',yerr=binned_shear_group_m_errs[m_i][shear_good] )
                ax.errorbar( binned_zs[m_i][magf_good], binned_magf_group_m_bests[m_i][magf_good],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r', markersize=3,
                             yerr=binned_magf_group_m_errs[m_i][magf_good] )
                ax.errorbar( binned_zs[m_i][overall_good]+z_shift, binned_overall_group_m_bests[m_i][overall_good],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*',yerr=binned_overall_group_m_errs[m_i][overall_good] )
                ax.set_ylim(13.1,15.9)
                ax.set_ylabel(r"$\log_{10}(M_{\rm gr}/M_{\rm sun})$",labelpad=5)
            elif(p_i==2):
                # Plot sat_frac
                ax = fig.add_subplot( num_m_bins, num_p_bins, p_i + num_p_bins*m_i + 1)
                ax.errorbar( binned_zs[m_i][shear_good]-z_shift, binned_shear_sat_frac_bests[m_i][shear_good],
                             color='b', linestyle='None', label="Shear fit",
                             marker='.',yerr=binned_shear_sat_frac_errs[m_i][shear_good] )
                ax.errorbar( binned_zs[m_i][magf_good], binned_magf_sat_frac_bests[m_i][magf_good],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r', markersize=3,
                             yerr=binned_magf_sat_frac_errs[m_i][magf_good] )
                ax.errorbar( binned_zs[m_i][overall_good]+z_shift, binned_overall_sat_frac_bests[m_i][overall_good],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*',yerr=binned_overall_sat_frac_errs[m_i][overall_good] )
                ax.set_ylim(0.01,0.49)
                ax.set_ylabel(r"$f_{\rm sat}$",labelpad=5)
            elif(p_i==3):
                # Plot sat_frac
                ax = fig.add_subplot( num_m_bins, num_p_bins, p_i + num_p_bins*m_i + 1)
                ax.errorbar( binned_zs[m_i][magf_good], binned_magf_Sigma_offset_bests[m_i][magf_good]/binned_Sigma_crits[m_i][magf_good],
                             color='r', linestyle='None', label="Magnification fit",
                             marker='o', markerfacecolor='none', markeredgecolor='r', markersize=3,
                             yerr=binned_magf_Sigma_offset_errs[m_i][magf_good]/binned_Sigma_crits[m_i][magf_good] )
                ax.errorbar( binned_zs[m_i][overall_good]+z_shift, binned_overall_Sigma_offset_bests[m_i][overall_good]/binned_Sigma_crits[m_i][overall_good],
                             color='k', linestyle='None', label="Overall fit",
                             marker='*',yerr=binned_overall_Sigma_offset_errs[m_i][overall_good]/binned_Sigma_crits[m_i][overall_good] )
                ax.set_ylim(-0.00549,0.00549)
                ax.set_ylabel(r"$\kappa_{\rm offset}$",labelpad=5)
                
            ax.set_xlim(0.2,1.3)
            
            # Label the mass
            xmin = 0.
            xmax = 1.
            ymin = 0.
            ymax = 1.
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.9, r"$M_{\rm mid}$=" + "%.1e" % m, size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            
            # set the labels and legend as appropriate
                            
            if(m_i!=num_m_bins-1): # Not on the bottom row
                ax.set_xticklabels([])
                
            if(m_i==num_m_bins-1): # bottom row
                ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
                ax.set_xticklabels([.2, .4, .6, .8, 1.0, 1.2],fontsize=12)
                
                if(p_i==0): # On bottom-left
                    ax.legend(loc='lower left', prop={'size':8}, numpoints=1)
    
    
    # Save the figure
    outfile_name = os.path.splitext(fitting_results_table_name)[0] + ".eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)