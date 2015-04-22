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

def main(argv):
    
    # Magic values
    
    high_z = False
    
    figsize = (8,4)
    labelsize = 8
    
    chi_2_i_min = 2
    
    default_lensing_signal_table = "/home/brg/git/CFHTLenS_cat/Data/gg_lensing_signal_20_bins_with_bf_models.dat"
    
    default_R_col_name = "shear_R_mean"
    
    default_z_col_name = "shear_lens_z_mean"
    if(not high_z):
        default_z_bins_min = 0.2
        default_z_bins_max = 0.7
        default_num_z_bins = 5
        z_str = ""
        xl_fontsize = 10
    else:
        default_z_bins_min = 0.7
        default_z_bins_max = 1.3
        default_num_z_bins = 6
        z_str = "_high_z"
        xl_fontsize = 8
    default_z_bins_log = False
    
    default_m_col_name = "shear_lens_m_mean"
    default_m_bins_min = 1e9
    default_m_bins_max = 1e12
    default_num_m_bins = 3
    default_m_bins_log = True
    
    default_mag_col_name = "shear_lens_mag_mean"
    default_mag_bins_min = 15
    default_mag_bins_max = 25
    default_num_mag_bins = 1
    default_mag_bins_log = False
    
    default_dS_col_name = "dS_t_mean"
    default_dS_err_col_name = "dS_t_stderr"
    
    default_Sigma_col_name = "Sigma"
    default_Sigma_err_col_name = "Sigma_stderr"
    
    default_bf_shear_model_dS_col_name = "bf_shear_model_dS_t"
    default_bf_shear_model_Sigma_col_name = "bf_shear_model_Sigma"
    
    default_bf_magf_model_dS_col_name = "bf_magf_model_dS_t"
    default_bf_magf_model_Sigma_col_name = "bf_magf_model_Sigma"
    
    default_bf_overall_model_dS_col_name = "bf_overall_model_dS_t"
    default_bf_overall_model_Sigma_col_name = "bf_overall_model_Sigma"
    
    default_shear_sigma_crit_col_name = "shear_Sigma_crit"
    default_magf_sigma_crit_col_name = "magf_Sigma_crit"
    
    default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    # Load in values from command line if provided, otherwise use defaults
    
    num_args = len(argv)
    cur_arg = 2
    
    if(num_args) < cur_arg:
        lensing_signal_table_name = default_lensing_signal_table
    else:
        lensing_signal_table_name = argv[cur_arg]
    
    # Load in the redshift table
    try:
        lensing_signal_table = ascii.read(lensing_signal_table_name)
    except:
        print("ERROR: Corr funcs table " + lensing_signal_table + " cannot be read.")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        R_col_name = default_R_col_name
    else:
        R_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        Rs = lensing_signal_table[R_col_name]
    except:
        print("ERROR: Cannot read column " + R_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_col_name = default_z_col_name
    else:
        z_col_name = argv[cur_arg]
        
    # Try to load in the redshift column
    try:
        redshifts = lensing_signal_table[z_col_name]
    except:
        print("ERROR: Cannot read column " + z_col_name + " from table " + lensing_signal_table_name + ".")
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
        masses = lensing_signal_table[m_col_name]
    except:
        print("ERROR: Cannot read column " + m_col_name + " from table " + lensing_signal_table_name + ".")
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
        _mags = lensing_signal_table[mag_col_name]
    except:
        print("ERROR: Cannot read column " + mag_col_name + " from table " + lensing_signal_table_name + ".")
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
        dS_col_name = default_dS_col_name
    else:
        dS_col_name = argv[cur_arg]
        
    # Try to load in the dS column
    try:
        dSs = lensing_signal_table[dS_col_name]
    except:
        print("ERROR: Cannot read column " + dS_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        dS_err_col_name = default_dS_err_col_name
    else:
        dS_err_col_name = argv[cur_arg]
        
    # Try to load in the dS_err column
    try:
        dS_errs = lensing_signal_table[dS_err_col_name]
    except:
        print("ERROR: Cannot read column " + dS_err_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        Sigma_col_name = default_Sigma_col_name
    else:
        Sigma_col_name = argv[cur_arg]
        
    # Try to load in the dS column
    try:
        Sigmas = lensing_signal_table[Sigma_col_name]
    except:
        print("ERROR: Cannot read column " + Sigma_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        Sigma_err_col_name = default_Sigma_err_col_name
    else:
        Sigma_err_col_name = argv[cur_arg]
        
    # Try to load in the dS_err column
    try:
        Sigma_errs = lensing_signal_table[Sigma_err_col_name]
    except:
        print("ERROR: Cannot read column " + Sigma_err_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_shear_model_dS_col_name = default_bf_shear_model_dS_col_name
    else:
        bf_shear_model_dS_col_name = argv[cur_arg]
        
    # Try to load in the bf_shear_model_dS column
    try:
        bf_shear_model_dSs = lensing_signal_table[bf_shear_model_dS_col_name]
    except:
        print("ERROR: Cannot read column " + bf_shear_model_dS_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_shear_model_Sigma_col_name = default_bf_shear_model_Sigma_col_name
    else:
        bf_shear_model_Sigma_col_name = argv[cur_arg]
        
    # Try to load in the bf_shear_model_Sigma column
    try:
        bf_shear_model_Sigmas = lensing_signal_table[bf_shear_model_Sigma_col_name]
    except:
        print("ERROR: Cannot read column " + bf_shear_model_Sigma_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_magf_model_dS_col_name = default_bf_magf_model_dS_col_name
    else:
        bf_magf_model_dS_col_name = argv[cur_arg]
        
    # Try to load in the bf_magf_model_dS column
    try:
        bf_magf_model_dSs = lensing_signal_table[bf_magf_model_dS_col_name]
    except:
        print("ERROR: Cannot read column " + bf_magf_model_dS_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_magf_model_Sigma_col_name = default_bf_magf_model_Sigma_col_name
    else:
        bf_magf_model_Sigma_col_name = argv[cur_arg]
        
    # Try to load in the bf_magf_model_Sigma column
    try:
        bf_magf_model_Sigmas = lensing_signal_table[bf_magf_model_Sigma_col_name]
    except:
        print("ERROR: Cannot read column " + bf_magf_model_Sigma_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_overall_model_dS_col_name = default_bf_overall_model_dS_col_name
    else:
        bf_overall_model_dS_col_name = argv[cur_arg]
        
    # Try to load in the bf_overall_model_dS column
    try:
        bf_overall_model_dSs = lensing_signal_table[bf_overall_model_dS_col_name]
    except:
        print("ERROR: Cannot read column " + bf_overall_model_dS_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_overall_model_Sigma_col_name = default_bf_overall_model_Sigma_col_name
    else:
        bf_overall_model_Sigma_col_name = argv[cur_arg]
        
    # Try to load in the bf_overall_model_Sigma column
    try:
        bf_overall_model_Sigmas = lensing_signal_table[bf_overall_model_Sigma_col_name]
    except:
        print("ERROR: Cannot read column " + bf_overall_model_Sigma_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        shear_sigma_crit_col_name = default_shear_sigma_crit_col_name
    else:
        shear_sigma_crit_col_name = argv[cur_arg]
        
    # Try to load in the bf_overall_model_Sigma column
    try:
        shear_sigma_crits = lensing_signal_table[shear_sigma_crit_col_name]
    except:
        print("ERROR: Cannot read column " + shear_sigma_crit_col_name + " from table " + lensing_signal_table_name + ".")
        return
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        magf_sigma_crit_col_name = default_magf_sigma_crit_col_name
    else:
        magf_sigma_crit_col_name = argv[cur_arg]
        
    # Try to load in the bf_overall_model_Sigma column
    try:
        magf_sigma_crits = lensing_signal_table[magf_sigma_crit_col_name]
    except:
        print("ERROR: Cannot read column " + magf_sigma_crit_col_name + " from table " + lensing_signal_table_name + ".")
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
    binned_dSs = []
    binned_dS_errs = []
    binned_Sigmas = []
    binned_neg_Sigmas = []
    binned_Sigma_errs = []
    binned_bf_shear_model_dSs = []
    binned_bf_shear_model_Sigmas = []
    binned_neg_bf_shear_model_Sigmas = []
    binned_bf_magf_model_dSs = []
    binned_bf_magf_model_Sigmas = []
    binned_neg_bf_magf_model_Sigmas = []
    binned_bf_overall_model_dSs = []
    binned_bf_overall_model_Sigmas = []
    binned_neg_bf_overall_model_Sigmas = []
    binned_shear_sigma_crits = []
    binned_magf_sigma_crits = []
    
    for z_i in xrange(num_z_bins):
        zR_list = []
        zdS_list = []
        zdS_err_list = []
        zSigma_list = []
        zneg_Sigma_list = []
        zSigma_err_list = []
        zbf_shear_model_dS_list = []
        zbf_shear_model_Sigma_list = []
        zneg_bf_shear_model_Sigma_list = []
        zbf_magf_model_dS_list = []
        zbf_magf_model_Sigma_list = []
        zneg_bf_magf_model_Sigma_list = []
        zbf_overall_model_dS_list = []
        zbf_overall_model_Sigma_list = []
        zneg_bf_overall_model_Sigma_list = []
        zshear_sigma_crits = []
        zmagf_sigma_crits = []
        for m_i in xrange(num_m_bins):
            zR_list.append([])
            zdS_list.append([])
            zdS_err_list.append([])
            zSigma_list.append([])
            zneg_Sigma_list.append([])
            zSigma_err_list.append([])
            zbf_shear_model_dS_list.append([])
            zbf_shear_model_Sigma_list.append([])
            zneg_bf_shear_model_Sigma_list.append([])
            zbf_magf_model_dS_list.append([])
            zbf_magf_model_Sigma_list.append([])
            zneg_bf_magf_model_Sigma_list.append([])
            zbf_overall_model_dS_list.append([])
            zbf_overall_model_Sigma_list.append([])
            zneg_bf_overall_model_Sigma_list.append([])
            zshear_sigma_crits.append([])
            zmagf_sigma_crits.append([])
        binned_Rs.append(zR_list)
        binned_dSs.append(zdS_list)
        binned_dS_errs.append(zdS_err_list)
        binned_Sigmas.append(zSigma_list)
        binned_neg_Sigmas.append(zneg_Sigma_list)
        binned_Sigma_errs.append(zSigma_err_list)
        binned_bf_shear_model_dSs.append(zbf_shear_model_dS_list)
        binned_bf_shear_model_Sigmas.append(zbf_shear_model_Sigma_list)
        binned_neg_bf_shear_model_Sigmas.append(zneg_bf_shear_model_Sigma_list)
        binned_bf_magf_model_dSs.append(zbf_magf_model_dS_list)
        binned_bf_magf_model_Sigmas.append(zbf_magf_model_Sigma_list)
        binned_neg_bf_magf_model_Sigmas.append(zneg_bf_magf_model_Sigma_list)
        binned_bf_overall_model_dSs.append(zbf_overall_model_dS_list)
        binned_bf_overall_model_Sigmas.append(zbf_overall_model_Sigma_list)
        binned_neg_bf_overall_model_Sigmas.append(zneg_bf_overall_model_Sigma_list)
        binned_shear_sigma_crits.append(zshear_sigma_crits)
        binned_magf_sigma_crits.append(zmagf_sigma_crits)
    
    for z, m, R, dS, dS_err, Sigma, Sigma_err, bf_shear_model_dS, bf_shear_model_Sigma, \
        bf_magf_model_dS, bf_magf_model_Sigma, bf_overall_model_dS, bf_overall_model_Sigma, \
        shear_sigma_crit, magf_sigma_crit in \
            zip(redshifts, masses, Rs, dSs, dS_errs, Sigmas, Sigma_errs, bf_shear_model_dSs, bf_shear_model_Sigmas, \
            bf_magf_model_dSs, bf_magf_model_Sigmas,bf_overall_model_dSs, bf_overall_model_Sigmas, \
            shear_sigma_crits, magf_sigma_crits):
        
        z_i = bf.get_bin_index(z, z_bins)
        if(z_i<0): continue
        
        m_i = bf.get_bin_index(m, m_bins)
        if(m_i<0): continue
        
        binned_Rs[z_i][m_i].append(R)
        binned_dSs[z_i][m_i].append(dS)
        binned_dS_errs[z_i][m_i].append(dS_err)
        binned_Sigmas[z_i][m_i].append(Sigma)
        binned_neg_Sigmas[z_i][m_i].append(-Sigma)
        binned_Sigma_errs[z_i][m_i].append(Sigma_err)
        binned_bf_shear_model_dSs[z_i][m_i].append(bf_shear_model_dS)
        binned_bf_shear_model_Sigmas[z_i][m_i].append(bf_shear_model_Sigma)
        binned_neg_bf_shear_model_Sigmas[z_i][m_i].append(-bf_shear_model_Sigma)
        binned_bf_magf_model_dSs[z_i][m_i].append(bf_magf_model_dS)
        binned_bf_magf_model_Sigmas[z_i][m_i].append(bf_magf_model_Sigma)
        binned_neg_bf_magf_model_Sigmas[z_i][m_i].append(-bf_magf_model_Sigma)
        binned_bf_overall_model_dSs[z_i][m_i].append(bf_overall_model_dS)
        binned_bf_overall_model_Sigmas[z_i][m_i].append(bf_overall_model_Sigma)
        binned_neg_bf_overall_model_Sigmas[z_i][m_i].append(-bf_overall_model_Sigma)
        binned_shear_sigma_crits[z_i][m_i].append(shear_sigma_crit)
        binned_magf_sigma_crits[z_i][m_i].append(magf_sigma_crit)
        
    # Convert all lists to arrays
    
    for z_i in xrange(num_z_bins):        
        for m_i in xrange(num_m_bins):
            binned_Rs[z_i][m_i] = np.array(binned_Rs[z_i][m_i])
            binned_dSs[z_i][m_i] = np.array(binned_dSs[z_i][m_i])
            binned_dS_errs[z_i][m_i] = np.array(binned_dS_errs[z_i][m_i])
            binned_Sigmas[z_i][m_i] = np.array(binned_Sigmas[z_i][m_i])
            binned_neg_Sigmas[z_i][m_i] = np.array(binned_neg_Sigmas[z_i][m_i])
            binned_Sigma_errs[z_i][m_i] = np.array(binned_Sigma_errs[z_i][m_i])
            binned_bf_shear_model_dSs[z_i][m_i] = np.array(binned_bf_shear_model_dSs[z_i][m_i])
            binned_bf_shear_model_Sigmas[z_i][m_i] = np.array(binned_bf_shear_model_Sigmas[z_i][m_i])
            binned_neg_bf_shear_model_Sigmas[z_i][m_i] = np.array(binned_neg_bf_shear_model_Sigmas[z_i][m_i])
            binned_bf_magf_model_dSs[z_i][m_i] = np.array(binned_bf_magf_model_dSs[z_i][m_i])
            binned_bf_magf_model_Sigmas[z_i][m_i] = np.array(binned_bf_magf_model_Sigmas[z_i][m_i])
            binned_neg_bf_magf_model_Sigmas[z_i][m_i] = np.array(binned_neg_bf_magf_model_Sigmas[z_i][m_i])
            binned_bf_overall_model_dSs[z_i][m_i] = np.array(binned_bf_overall_model_dSs[z_i][m_i])
            binned_bf_overall_model_Sigmas[z_i][m_i] = np.array(binned_bf_overall_model_Sigmas[z_i][m_i])
            binned_neg_bf_overall_model_Sigmas[z_i][m_i] = np.array(binned_neg_bf_overall_model_Sigmas[z_i][m_i])
            binned_shear_sigma_crits[z_i][m_i] = np.array(binned_shear_sigma_crits[z_i][m_i])
            binned_magf_sigma_crits[z_i][m_i] = np.array(binned_magf_sigma_crits[z_i][m_i])
            
    # Set up chi-squared storage arrays
    shear_dS_chi2s = np.empty((num_z_bins,num_m_bins))
    magf_dS_chi2s = np.empty((num_z_bins,num_m_bins))
    shear_Sigma_chi2s = np.empty((num_z_bins,num_m_bins))
    magf_Sigma_chi2s = np.empty((num_z_bins,num_m_bins))
        
    # Do the shear plot now
    fig = pyplot.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("Projected Separation (kpc)",labelpad=10)
    ax.set_ylabel(r"$\gamma_{\rm t}$",labelpad=25)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    for z_i in xrange(num_z_bins):
        
        z = z_bins_mids[z_i]
        
        for m_i in xrange(num_m_bins):
        
            m = m_bins_mids[m_i]
            
            # Calculate the chi-squared values for this bin            
            shear_chi2s = np.square((binned_dSs[z_i][m_i] - binned_bf_shear_model_dSs[z_i][m_i])/binned_dS_errs[z_i][m_i])
            magf_chi2s = np.square((binned_dSs[z_i][m_i] - binned_bf_magf_model_dSs[z_i][m_i])/binned_dS_errs[z_i][m_i])
            
            shear_dS_chi2s[z_i,m_i] = np.sum(shear_chi2s[chi_2_i_min:])
            magf_dS_chi2s[z_i,m_i] = np.sum(magf_chi2s[chi_2_i_min:])
            
            # Plot this bin
            ax = fig.add_subplot( num_m_bins, num_z_bins, z_i + num_z_bins*m_i + 1)
            ax.set_yscale("log", nonposy='clip')
            ax.errorbar( binned_Rs[z_i][m_i], binned_dSs[z_i][m_i]/binned_shear_sigma_crits[z_i][m_i],
                         color='b', linestyle='None', label="Measured values",
                         marker='.',yerr=binned_dS_errs[z_i][m_i]/binned_shear_sigma_crits[z_i][m_i] )
            ax.plot( binned_Rs[z_i][m_i], binned_bf_shear_model_dSs[z_i][m_i]/binned_shear_sigma_crits[z_i][m_i],
                     'k', linestyle='--', linewidth=2,
                     label="Best fit to shear only")
            ax.plot( binned_Rs[z_i][m_i], binned_bf_magf_model_dSs[z_i][m_i]/binned_shear_sigma_crits[z_i][m_i],
                     'r', linestyle=':', linewidth=2,
                     label="Best fit to mag. only")
            #ax.legend( [emptiness,emptiness] , [r"$z_{mid}$="+ str(z) ,r"$M_{mid}$=" + "%.1E" % m],loc='upper right')
            ax.set_ylim(0.0003,0.03)
            
            # Label the redshift and mass
            xmin = 0.
            xmax = 1.
            ymin = 0.
            ymax = 1.
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.9, r"$z_{mid}$="+ str(z), size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.8, r"$M_{mid}$=" + "%.1e" % m, size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.7, r"$\chi^2$=" +
                     "%2.1f, %2.1f" % (shear_dS_chi2s[z_i,m_i], magf_dS_chi2s[z_i,m_i]), size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            
            # set the labels as appropriate
            if(z_i!=0): # Not on the left column
                ax.set_yticklabels([])
                
            if(m_i!=num_m_bins-1): # Not on the bottom row
                ax.set_xticklabels([])
                
            if((z_i==0) and (m_i==num_m_bins-1)): # bottom-left
                ax.set_yticks([0.001, 0.01])
                ax.set_yticklabels([0.001, 0.01],fontsize=10)
                ax.set_xticks([0, 500,1000,1500,2000])
                ax.set_xticklabels([0, 500,1000,1500,2000],fontsize=xl_fontsize)
                continue
                
            if(z_i==0): # left column
                ax.set_yticks([0.001, 0.01])
                ax.set_yticklabels([0.001, 0.01],fontsize=10)
                
            if(m_i==num_m_bins-1): # bottom row
                ax.set_xticks([500,1000,1500,2000])
                ax.set_xticklabels([500,1000,1500,2000],fontsize=xl_fontsize)
    
    
    # Save the figure
    outfile_name = os.path.splitext(lensing_signal_table_name)[0] + z_str + "_gamma.eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()
    
    # Print stats on the chi-squared value
    print("shear dS Chi^2 = " + str(np.mean(shear_dS_chi2s)) + " +/- " +
          str(np.std(shear_dS_chi2s)/np.sqrt(np.size(shear_dS_chi2s)-1)) )
    print("magf dS Chi^2 = " + str(np.mean(magf_dS_chi2s)) + " +/- " +
          str(np.std(magf_dS_chi2s)/np.sqrt(np.size(magf_dS_chi2s)-1)) )
    print("")
        
    # Do the mag plot now
    fig = pyplot.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("Projected Separation (kpc)",labelpad=20)
    ax.set_ylabel(r"$\left|\kappa\right|$",labelpad=25)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    for z_i in xrange(num_z_bins):
        
        z = z_bins_mids[z_i]
        
        for m_i in xrange(num_m_bins):
        
            m = m_bins_mids[m_i]
            
            # Calculate the chi-squared values for this bin            
            shear_chi2s = np.square((binned_Sigmas[z_i][m_i] - binned_bf_shear_model_Sigmas[z_i][m_i])/binned_Sigma_errs[z_i][m_i])
            magf_chi2s = np.square((binned_Sigmas[z_i][m_i] - binned_bf_magf_model_Sigmas[z_i][m_i])/binned_Sigma_errs[z_i][m_i])
            
            shear_Sigma_chi2s[z_i,m_i] = np.sum(shear_chi2s[chi_2_i_min:])
            magf_Sigma_chi2s[z_i,m_i] = np.sum(magf_chi2s[chi_2_i_min:])
    
            ax = fig.add_subplot( num_m_bins, num_z_bins, z_i + num_z_bins*m_i + 1)
            ax.set_yscale("log", nonposy='clip')
            ax.errorbar( binned_Rs[z_i][m_i], binned_Sigmas[z_i][m_i]/binned_magf_sigma_crits[z_i][m_i],
                         color='r', linestyle='None', label="Measured values",
                         marker='.', yerr=binned_Sigma_errs[z_i][m_i]/binned_shear_sigma_crits[z_i][m_i] )
            ax.errorbar( binned_Rs[z_i][m_i], binned_neg_Sigmas[z_i][m_i]/binned_magf_sigma_crits[z_i][m_i],
                         color='r', linestyle='None', label="Negative Measured values",
                         marker='o', markerfacecolor='none', markeredgecolor='r',
                         yerr=binned_Sigma_errs[z_i][m_i]/binned_shear_sigma_crits[z_i][m_i] )
            ax.plot( binned_Rs[z_i][m_i], binned_bf_shear_model_Sigmas[z_i][m_i]/binned_magf_sigma_crits[z_i][m_i],
                     'b', linestyle='--',  linewidth=2,
                     label="Best fit to shear only")
            ax.plot( binned_Rs[z_i][m_i], binned_neg_bf_shear_model_Sigmas[z_i][m_i]/binned_magf_sigma_crits[z_i][m_i],
                     'b', linestyle='--', linewidth=2,
                     label="Best fit to shear only")
            ax.plot( binned_Rs[z_i][m_i], binned_bf_magf_model_Sigmas[z_i][m_i]/binned_magf_sigma_crits[z_i][m_i],
                     'k', linestyle=':', linewidth=2,
                     label="Best fit to mag. only")
            ax.plot( binned_Rs[z_i][m_i], binned_neg_bf_magf_model_Sigmas[z_i][m_i]/binned_magf_sigma_crits[z_i][m_i],
                     'k', linestyle=':', linewidth=2,
                     label="Best fit to mag. only")
            #ax.legend( [emptiness,emptiness] , [r"$z_{mid}$="+ str(z) ,r"$M_{mid}$=" + "%.1E" % m],loc='upper right')
            ax.set_ylim(0.0003,0.03)
            
            # Label the redshift and mass
            xmin = 0.
            xmax = 1.
            ymin = 0.
            ymax = 1.
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.9, r"$z_{mid}$="+ str(z), size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.8, r"$M_{mid}$=" + "%.1E" % m, size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            ax.text(xmin+(xmax-xmin)*0.95, ymin+(ymax-ymin)*0.7, r"$\chi^2$=" +
                    "%2.1f, %2.1f" % (magf_Sigma_chi2s[z_i,m_i], shear_Sigma_chi2s[z_i,m_i]), size=labelsize,
                    horizontalalignment='right', transform = ax.transAxes)
            
            # set the labels as appropriate
            if(z_i!=0): # Not on the left column
                ax.set_yticklabels([])
                
            if(m_i!=num_m_bins-1): # Not on the bottom row
                ax.set_xticklabels([])
                
            if((z_i==0) and (m_i==num_m_bins-1)): # bottom-left
                ax.set_yticks([0.001, 0.01])
                ax.set_yticklabels([0.001, 0.01],fontsize=10)
                ax.set_xticks([0, 500,1000,1500,2000])
                ax.set_xticklabels([0, 500,1000,1500,2000],fontsize=xl_fontsize)
                continue
                
            if(z_i==0): # left column
                ax.set_yticks([0.001, 0.01])
                ax.set_yticklabels([0.001, 0.01],fontsize=10)
                
            if(m_i==num_m_bins-1): # bottom row
                ax.set_xticks([500,1000,1500,2000])
                ax.set_xticklabels([500,1000,1500,2000],fontsize=xl_fontsize)
    
    
    # Save the figure
    outfile_name = os.path.splitext(lensing_signal_table_name)[0] + z_str + "_kappa.eps"
    pyplot.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
        
    # Show the figure
    pyplot.show()
    
    # Print stats on the chi-squared value
    print("magf Sigma Chi^2 = " + str(np.mean(magf_Sigma_chi2s)) + " +/- " +
          str(np.std(magf_Sigma_chi2s)/np.sqrt(np.size(magf_Sigma_chi2s)-1)) )
    print("shear Sigma Chi^2 = " + str(np.mean(shear_Sigma_chi2s)) + " +/- " +
          str(np.std(shear_Sigma_chi2s)/np.sqrt(np.size(shear_Sigma_chi2s)-1)) )

if __name__ == "__main__":
    main(sys.argv)