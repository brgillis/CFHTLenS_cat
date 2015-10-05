#!/usr/bin/env python

""" plot_covariance.py

    Created 23 Sep 2015

    Constructs a covariance triangle plot from shear/magnification fitting
    data.

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess as sbp
import sys

import matplotlib

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

# Magic values
bins_to_plot = ((0.2,10), (1,10))

paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
data_name_root = "/home/brg/git/CFHTLenS_cat/Data/gg_lensing_signal_20_bins_fitting_results"
    
# Plottable ranges for each parameter
sat_m_range = (10.5,13.499)
group_m_range = (11.,15.99)
sat_frac_range = (0.,0.99)
kappa_offset_range = (-4,1.99)    

# Labels for each parameter
sat_m_label = r"$\log_{\rm 10}\left( M_{\rm 1h}/M_{\rm sun} \right)$"
group_m_label = r"$\log_{\rm 10}\left( M_{\rm gr}/M_{\rm sun} \right)$"
sat_frac_label = r"$f_{\rm sat}$"
kappa_offset_label = r"$\kappa_{\rm offset}/10^{-3}$"

def main(argv):
    
    for bin_to_plot in bins_to_plot:
        z_tag = str(bin_to_plot[0]).replace(".","p")
        m_tag = str(bin_to_plot[1]).replace(".","p")
        
        data_name_tag = "_z-" + z_tag + "_m-" + m_tag + "_"
    
        shear_data = np.genfromtxt(data_name_root + data_name_tag + 'shear_test_points.dat',
                             names='sat_m,group_m,sat_frac,kappa_offset',
                             usecols=(0,1,2,3))
        magf_data = np.genfromtxt(data_name_root + data_name_tag + 'magf_test_points.dat',
                             names='sat_m,group_m,sat_frac,kappa_offset',
                             usecols=(0,1,2,3))
        overall_data = np.genfromtxt(data_name_root + data_name_tag + 'overall_test_points.dat',
                             names='sat_m,group_m,sat_frac,kappa_offset',
                             usecols=(0,1,2,3))
        
        sat_m = shear_data['sat_m'], magf_data['sat_m'], overall_data['sat_m']
        group_m = shear_data['group_m'], magf_data['group_m'], overall_data['group_m']
        sat_frac = shear_data['sat_frac'], magf_data['sat_frac'], overall_data['sat_frac']
        kappa_offset = shear_data['kappa_offset']*1.0e3, magf_data['kappa_offset']*1.0e3, \
            overall_data['kappa_offset']*1.0e3
        
        # Number of bins in 2d plots (per axis)
        num_bins_2d = 20
        
        # Number of bins in 1d plots
        num_bins_1d = 50
        
        fig = plt.figure(4, (8,8))
        fig.subplots_adjust(hspace=0.001, wspace=0.001, left=0.08, bottom=0.095, top=0.975, right=0.93)
        
        # Set up each of the 2D plots
        for  x_data,   x_range,        x_label,        y_data,       y_range, \
                y_label,     plot_shear, pos, bottom_row, left_col in (
            (sat_m,    sat_m_range,    sat_m_label,    group_m,      group_m_range,
                group_m_label,      True,  13, True,  True),
            (sat_m,    sat_m_range,    sat_m_label,    sat_frac,     sat_frac_range,
                sat_frac_label,     True,  9, False,  True),
            (sat_m,    sat_m_range,    sat_m_label,    kappa_offset, kappa_offset_range,
                kappa_offset_label, False, 5, False,  True),
            (kappa_offset,  kappa_offset_range,  kappa_offset_label,  group_m,     group_m_range,
                group_m_label,     False,  14,  True, False),
            (kappa_offset,  kappa_offset_range,  kappa_offset_label,  sat_frac, sat_frac_range,
                sat_frac_label, False, 10,  False, False),
            (sat_frac, sat_frac_range, sat_frac_label, group_m, group_m_range,
                group_m_label, True, 15,  True, False)):
            
            # Work on the plot in this position
            ax = plt.subplot(4,4,pos)
            
            # Convert to 2d histogram for each of the shear, magnification, and overall data
            for index, label, color, linestyle, lineweight in ((0, "Shear",         "b", "--", 1),
                                                   (1, "Magnification", "r", "-", 1),
                                                   (2, "Combined",      "k", "-", 1)):
                # Skip shear if necessary
                if((not plot_shear) and (index==0)):
                    continue
                
                hist2D, xedges, yedges = np.histogram2d(x_data[index], y_data[index],
                                                        bins=[num_bins_2d,num_bins_2d],
                                                        range=[x_range,y_range],
                                                        normed=False)
            
                # Plot Monte-Carlo samples as 2D histogram.
                hist2D = np.transpose(hist2D)  # Beware: np switches axes, so switch back.
            
                extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            
                # We have to find the proper levels for 68% and 95% of the data
                totcount = hist2D.sum()
                levelstep = 0.0001*np.max(hist2D)
                
                #95% first
                goalcount = 0.95*totcount
                curcount = totcount
                
                testlevel = levelstep
                while(curcount>goalcount):
                    curcount = hist2D[hist2D > testlevel].sum()
                    testlevel += levelstep
                
                level95 = testlevel - levelstep/2
                
                #68% now
                goalcount = 0.68*totcount
                curcount = totcount
                
                while(curcount>goalcount):
                    curcount = hist2D[hist2D > testlevel].sum()
                    testlevel += levelstep
                
                level68 = testlevel - levelstep/2
                
                levels = (level68, level95)
                cset = plt.contour(hist2D, levels, origin='lower',colors=[color,color],
                                   linewidths=(2*lineweight, lineweight),extent=extent,label=label)
                for c in cset.collections:
                    c.set_linestyle(linestyle)
                
                # Add to legend only for the plot in bottom-left (so we don't get duplicates)
                if(bottom_row and left_col):
                    _proxy = [plt.plot(0,0,'-',lw=3, color=color,linestyle=linestyle,label=label)]
        
                    # Plot the legend for these plots in the upper right corner
                    plt.legend( loc='upper right',bbox_to_anchor=(-0.2, -0.3, 1, 1),
                                bbox_transform=plt.gcf().transFigure)
            
            plt.xlim(x_range)
            plt.ylim(y_range)
            
            # If on the bottom row, set the x label
            if(bottom_row):
                plt.xlabel(x_label,fontsize=14)
            else:
                # Else no x ticks
                ax.set_xticklabels([])
                
            # If on the left edge, set the y label
            if(left_col):
                plt.ylabel(y_label,fontsize=14)
            else:
                # Else no y tick labels
                ax.set_yticklabels([])
                
            # Set the font size for the ticks
            plt.yticks(fontsize=10)
            plt.xticks(fontsize=10)
        
        # Now do the 1D histograms
        for  data,         data_range,         data_label,    plot_shear, pos, bottom_row, left_col in (
            (sat_m,        sat_m_range,        sat_m_label,   True,       1,  False,        True),
            (group_m,      group_m_range,      group_m_label, True,       16,  True,       False),
            (sat_frac,     sat_frac_range,     sat_frac_label,   True,       11,   False,       False),
            (kappa_offset, kappa_offset_range, kappa_offset_label,   False,      6,   False,       False)):
            
            # Set the position we're plotting now
            ax = plt.subplot(4,4,pos)
    
            # Determine the peak for all fitting types
            peak = 0.
    
            # Convert to 1d histogram for each of the shear, magnification, and overall data
            for index, label, color, linestyle, lineweight in ((0, "Shear",         "b", "--", 1),
                                                   (1, "Magnification", "r", "-", 1),
                                                   (2, "Combined",      "k", "-", 1)):
                # Skip shear if necessary
                if((not plot_shear) and (index==0)):
                    continue
    
                hist_height = np.histogram(data[index], bins=num_bins_1d, range=data_range, normed=True)[0]
                # Restore positions lost by binning.
                hist_data = data_range[0] + (data_range[1]-data_range[0]) * \
                    np.array(range(0,len(hist_height)))/float(len(hist_height)-1)
    
                # bottom right panel: projected density of x.
                plt.plot(hist_data, hist_height, color=color, linestyle=linestyle,
                         linewidth=lineweight)
                
                cur_peak = np.max(hist_height)
                if(cur_peak>peak):
                    peak = cur_peak
            
            plt.yticks([])
            plt.xlim(data_range)
            plt.ylim(0.0, 1.1*peak)
            
            # Add an x label if on the bottom row
            if(bottom_row):
                plt.xlabel(data_label,fontsize=14)
            else:
                # Else no x ticks
                ax.set_xticklabels([])
            
            # Add a y label if on the left column
            if(left_col):
                plt.ylabel(data_label,fontsize=14)
            else:
                # Else no y tick labels
                ax.set_yticklabels([])
                
            plt.xticks(fontsize=10)
            
        # Set up the label for the whole plot
        ax = fig.add_subplot(1,1,1)
        ax.spines["top"].set_color("none")
        ax.spines["bottom"].set_color("none")
        ax.spines["left"].set_color("none")
        ax.spines["right"].set_color("none")
        ax.tick_params(labelcolor="w", top="off", bottom="off", left="off", right="off")
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_axis_bgcolor("none")
        
        z_min = bin_to_plot[0]
        z_max = bin_to_plot[0] + 0.1
        
        m_min = bin_to_plot[1]
        m_max = bin_to_plot[1] + 1
        
        plot_label = r"$" + str(z_min) + r" \leq z_{\rm lens} < " + str(z_max) + "$, " + \
            r"$\;\;10^{" + str(m_min) + r"} M_{\rm sun} \leq M_{\rm lens} < 10^{" + \
            str(m_max) + r"} M_{\rm sun}$"
            
        ax.set_xlabel(plot_label,labelpad=30,fontsize=16)
    
        # Save the figure
        outfile_name = data_name_root + data_name_tag + "_covar.eps"
        plt.savefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
        
        # Copy it to the paper location
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
            
        # Show the figure
        plt.show()

if __name__ == "__main__":
    main(sys.argv)
