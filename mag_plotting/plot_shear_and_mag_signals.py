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


def main(argv):
    
    # Magic values
    
    high_z = True
    
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
    
    default_magf_sigma_offset_col_name = "bf_magf_Sigma_offset"
    default_overall_sigma_offset_col_name = "bf_overall_Sigma_offset"
    
    default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    # Load in values from command line if provided, otherwise use defaults
    
    reader = cr.columns_reader()
    
    num_args = len(argv)
    cur_arg = 2
    
    if(num_args) < cur_arg:
        lensing_signal_table_name = default_lensing_signal_table
    else:
        lensing_signal_table_name = argv[cur_arg]
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        R_col_name = default_R_col_name
    else:
        R_col_name = argv[cur_arg]
        
    reader.add("R",R_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        z_col_name = default_z_col_name
    else:
        z_col_name = argv[cur_arg]
        
    reader.add("z",z_col_name)
        
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
        
    reader.add("m",m_col_name)
        
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
        
    reader.add("mag",mag_col_name)
        
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
        
    reader.add("dS",dS_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        dS_err_col_name = default_dS_err_col_name
    else:
        dS_err_col_name = argv[cur_arg]
        
    reader.add("dS_err",dS_err_col_name)

        
    cur_arg += 1
    if(num_args) <= cur_arg:
        Sigma_col_name = default_Sigma_col_name
    else:
        Sigma_col_name = argv[cur_arg]
        
    reader.add("Sigma",Sigma_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        Sigma_err_col_name = default_Sigma_err_col_name
    else:
        Sigma_err_col_name = argv[cur_arg]
        
    reader.add("Sigma_err",Sigma_err_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_shear_model_dS_col_name = default_bf_shear_model_dS_col_name
    else:
        bf_shear_model_dS_col_name = argv[cur_arg]
        
    reader.add("bf_shear_model_dS",bf_shear_model_dS_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_shear_model_Sigma_col_name = default_bf_shear_model_Sigma_col_name
    else:
        bf_shear_model_Sigma_col_name = argv[cur_arg]
        
    reader.add("bf_shear_model_Sigma",bf_shear_model_Sigma_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_magf_model_dS_col_name = default_bf_magf_model_dS_col_name
    else:
        bf_magf_model_dS_col_name = argv[cur_arg]
        
    reader.add("bf_magf_model_dS",bf_magf_model_dS_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_magf_model_Sigma_col_name = default_bf_magf_model_Sigma_col_name
    else:
        bf_magf_model_Sigma_col_name = argv[cur_arg]
        
    reader.add("bf_magf_model_Sigma",bf_magf_model_Sigma_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_overall_model_dS_col_name = default_bf_overall_model_dS_col_name
    else:
        bf_overall_model_dS_col_name = argv[cur_arg]
        
    reader.add("bf_overall_model_dS",bf_overall_model_dS_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        bf_overall_model_Sigma_col_name = default_bf_overall_model_Sigma_col_name
    else:
        bf_overall_model_Sigma_col_name = argv[cur_arg]
        
    reader.add("bf_overall_model_Sigma",bf_overall_model_Sigma_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        shear_sigma_crit_col_name = default_shear_sigma_crit_col_name
    else:
        shear_sigma_crit_col_name = argv[cur_arg]
        
    reader.add("shear_sigma_crit",shear_sigma_crit_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        magf_sigma_crit_col_name = default_magf_sigma_crit_col_name
    else:
        magf_sigma_crit_col_name = argv[cur_arg]
        
    reader.add("magf_sigma_crit",magf_sigma_crit_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        magf_sigma_offset_col_name = default_magf_sigma_offset_col_name
    else:
        magf_sigma_offset_col_name = argv[cur_arg]
        
    reader.add("magf_sigma_offset",magf_sigma_offset_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        overall_sigma_offset_col_name = default_overall_sigma_offset_col_name
    else:
        overall_sigma_offset_col_name = argv[cur_arg]
        
    reader.add("overall_sigma_offset",overall_sigma_offset_col_name)
        
    cur_arg += 1
    if(num_args) <= cur_arg:
        paper_location = default_paper_location
    else:
        paper_location = argv[cur_arg]
        
    cols = reader.read(lensing_signal_table_name)
    
    binned_cols = []
           
    for _col_i in xrange(reader.num_cols()):
        binned_cols.append([])
    
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
    
    for col_i in xrange(reader.num_cols()):
        for _z_i in xrange(num_z_bins):
            new_list = []
            for _m_i in xrange(num_m_bins):
                new_list.append([])
            binned_cols[col_i].append(new_list)

    for elems in zip(*cols):
        
        z = elems[reader.index("z")]
        m = elems[reader.index("m")]
        
        z_i = bf.get_bin_index(z+0.001, z_bins)
        if(z_i<0): continue
        m_i = bf.get_bin_index(m+0.001, m_bins)
        if(m_i<0): continue
        
        for col_i in xrange(reader.num_cols()):
            binned_cols[col_i][z_i][m_i].append(elems[col_i])
        
    # Convert all lists to arrays
    for col_i in xrange(reader.num_cols()):
        for z_i in xrange(num_z_bins):
            for m_i in xrange(num_m_bins):
                binned_cols[col_i][z_i][m_i] = np.array(binned_cols[col_i][z_i][m_i])
            
    # Set up chi-squared storage arrays
    shear_dS_chi2s = np.empty((num_z_bins,num_m_bins))
    magf_dS_chi2s = np.empty((num_z_bins,num_m_bins))
    overall_dS_chi2s = np.empty((num_z_bins,num_m_bins))
    shear_Sigma_chi2s = np.empty((num_z_bins,num_m_bins))
    magf_Sigma_chi2s = np.empty((num_z_bins,num_m_bins))
    overall_Sigma_chi2s = np.empty((num_z_bins,num_m_bins))
        
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
            shear_chi2s = np.square((binned_cols[reader.index("dS")][z_i][m_i] - \
                                      binned_cols[reader.index("bf_shear_model_dS")][z_i][m_i])/ \
                                        binned_cols[reader.index("dS_err")][z_i][m_i])
            magf_chi2s = np.square((binned_cols[reader.index("dS")][z_i][m_i] - \
                                     binned_cols[reader.index("bf_magf_model_dS")][z_i][m_i])/ \
                                        binned_cols[reader.index("dS_err")][z_i][m_i])
            overall_chi2s = np.square((binned_cols[reader.index("dS")][z_i][m_i] - \
                                       binned_cols[reader.index("bf_overall_model_dS")][z_i][m_i])/ \
                                        binned_cols[reader.index("dS_err")][z_i][m_i])
            
            shear_dS_chi2s[z_i,m_i] = np.sum(shear_chi2s[chi_2_i_min:])
            magf_dS_chi2s[z_i,m_i] = np.sum(magf_chi2s[chi_2_i_min:])
            overall_dS_chi2s[z_i,m_i] = np.sum(overall_chi2s[chi_2_i_min:])
            
            # Plot this bin
            ax = fig.add_subplot( num_m_bins, num_z_bins, z_i + num_z_bins*m_i + 1)
            ax.set_yscale("log", nonposy='clip')
            ax.errorbar( binned_cols[reader.index("R")][z_i][m_i], binned_cols[reader.index("dS")][z_i][m_i]/binned_cols[reader.index("shear_sigma_crit")][z_i][m_i],
                         color='b', linestyle='None', label="Measured values",
                         marker='.',yerr=binned_cols[reader.index("dS_err")][z_i][m_i]/binned_cols[reader.index("shear_sigma_crit")][z_i][m_i] )
            ax.plot( binned_cols[reader.index("R")][z_i][m_i], binned_cols[reader.index("bf_shear_model_dS")][z_i][m_i]/binned_cols[reader.index("shear_sigma_crit")][z_i][m_i],
                     'b', linestyle='--', linewidth=2,
                     label="Best fit to shear only")
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
                 "%2.1f" % shear_dS_chi2s[z_i,m_i], size=labelsize,
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
    pyplot.writefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
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
    print("overall dS Chi^2 = " + str(np.mean(overall_dS_chi2s)) + " +/- " +
          str(np.std(overall_dS_chi2s)/np.sqrt(np.size(overall_dS_chi2s)-1)) )
    print("")
        
    # Do the mag plot now
    fig = pyplot.figure(figsize=figsize)
    fig.subplots_adjust(wspace=0, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("Projected Separation (kpc)",labelpad=10)
    ax.set_ylabel(r"$\left|\kappa\right|$",labelpad=25)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    for z_i in xrange(num_z_bins):
        
        z = z_bins_mids[z_i]
        
        for m_i in xrange(num_m_bins):
        
            m = m_bins_mids[m_i]
            
            # Calculate the chi-squared values for this bin            
            shear_chi2s = np.square((binned_cols[reader.index("Sigma")][z_i][m_i] - binned_cols[reader.index("bf_shear_model_Sigma")][z_i][m_i])/binned_cols[reader.index("Sigma_err")][z_i][m_i])
            magf_chi2s = np.square((binned_cols[reader.index("Sigma")][z_i][m_i] - binned_cols[reader.index("bf_magf_model_Sigma")][z_i][m_i])/binned_cols[reader.index("Sigma_err")][z_i][m_i])
            overall_chi2s = np.square((binned_cols[reader.index("Sigma")][z_i][m_i] - binned_cols[reader.index("bf_overall_model_Sigma")][z_i][m_i])/binned_cols[reader.index("Sigma_err")][z_i][m_i])
            
            shear_Sigma_chi2s[z_i,m_i] = np.sum(shear_chi2s[chi_2_i_min:])
            magf_Sigma_chi2s[z_i,m_i] = np.sum(magf_chi2s[chi_2_i_min:])
            overall_Sigma_chi2s[z_i,m_i] = np.sum(overall_chi2s[chi_2_i_min:])
    
            ax = fig.add_subplot( num_m_bins, num_z_bins, z_i + num_z_bins*m_i + 1)
            ax.set_yscale("symlog",linthreshy=0.001)
            ax.errorbar( binned_cols[reader.index("R")][z_i][m_i], binned_cols[reader.index("Sigma")][z_i][m_i]/binned_cols[reader.index("magf_sigma_crit")][z_i][m_i],
                         color='r', linestyle='None', label="Measured values",
                         marker='.', yerr=binned_cols[reader.index("Sigma_err")][z_i][m_i]/binned_cols[reader.index("shear_sigma_crit")][z_i][m_i] )
            ax.plot( binned_cols[reader.index("R")][z_i][m_i], binned_cols[reader.index("bf_shear_model_Sigma")][z_i][m_i]/binned_cols[reader.index("magf_sigma_crit")][z_i][m_i],
                     'b', linestyle='--',  linewidth=2,
                     label="Best fit to shear only")
            ax.plot( binned_cols[reader.index("R")][z_i][m_i], binned_cols[reader.index("bf_magf_model_Sigma")][z_i][m_i]/binned_cols[reader.index("magf_sigma_crit")][z_i][m_i],
                     'r', linestyle='-', linewidth=1,
                     label="Best fit to mag. only")
            ax.plot( binned_cols[reader.index("R")][z_i][m_i], binned_cols[reader.index("bf_overall_model_Sigma")][z_i][m_i]/binned_cols[reader.index("magf_sigma_crit")][z_i][m_i],
                     'k', linestyle='-.', linewidth=2,
                     label="Best fit to both")
            
            # Plot the zero line
            ax.plot([0.,2000.],[0.,0.],label=None,color="k",linestyle="dashed")
            
            ax.set_ylim(-0.011,0.09)
            
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
                     "%2.1f, %2.1f, %2.1f" % (shear_Sigma_chi2s[z_i,m_i],
                                              magf_Sigma_chi2s[z_i,m_i],
                                              overall_Sigma_chi2s[z_i,m_i]),
                    size=labelsize, horizontalalignment='right', transform = ax.transAxes)
            
            # set the labels as appropriate
            if(z_i!=0): # Not on the left column
                ax.set_yticklabels([])
                
            if(m_i!=num_m_bins-1): # Not on the bottom row
                ax.set_xticklabels([])
                
            yticks = [-0.01,-0.001,0.,0.001, 0.01]
            xticks = [0, 500,1000,1500,2000]
            xticks_br = [500,1000,1500,2000]
                
            if((z_i==0) and (m_i==num_m_bins-1)): # bottom-left
                ax.set_yticks(yticks)
                ax.set_yticklabels(yticks,fontsize=10)
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticks,fontsize=xl_fontsize)
                continue
                
            if(z_i==0): # left column
                ax.set_yticks(yticks)
                ax.set_yticklabels(yticks,fontsize=10)
                
            if(m_i==num_m_bins-1): # bottom row
                ax.set_xticks(xticks_br)
                ax.set_xticklabels(xticks_br,fontsize=xl_fontsize)
    
    
    # Save the figure
    outfile_name = os.path.splitext(lensing_signal_table_name)[0] + z_str + "_kappa.eps"
    pyplot.writefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
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
    print("overall Sigma Chi^2 = " + str(np.mean(overall_Sigma_chi2s)) + " +/- " +
          str(np.std(overall_Sigma_chi2s)/np.sqrt(np.size(overall_Sigma_chi2s)-1)) )

if __name__ == "__main__":
    main(sys.argv)