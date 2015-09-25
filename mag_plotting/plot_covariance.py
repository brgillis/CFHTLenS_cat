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
import matplotlib.gridspec as gridspec
import sys

from matplotlib.ticker import FixedLocator
from matplotlib import rc

rc('text', usetex=True)

# Magic values
bins_to_plot = ((0.2,10.), (0.9, 9.), (1.0,10.))

def main(argv):
    shear_data = np.genfromtxt('gg_lensing_signal_20_bins_fitting_results_z-0p2_m-10_shear_test_points.dat',
                         names='sat_m,group_m,sat_frac,Sigma_offset',
                         usecols=(0,1,2,3))
    magf_data = np.genfromtxt('gg_lensing_signal_20_bins_fitting_results_z-0p2_m-10_magf_test_points.dat',
                         names='sat_m,group_m,sat_frac,Sigma_offset',
                         usecols=(0,1,2,3))
    overall_data = np.genfromtxt('gg_lensing_signal_20_bins_fitting_results_z-0p2_m-10_overall_test_points.dat',
                         names='sat_m,group_m,sat_frac,Sigma_offset',
                         usecols=(0,1,2,3))
    
    sat_m = shear_data['sat_m'], magf_data['sat_m'], overall_data['sat_m']
    group_m = shear_data['group_m'], magf_data['group_m'], overall_data['group_m']
    sat_frac = shear_data['sat_frac'], magf_data['sat_frac'], overall_data['sat_frac']
    Sigma_offset = shear_data['Sigma_offset'], magf_data['Sigma_offset'], overall_data['Sigma_offset']
    
    # Plottable ranges for each parameter
    sat_m_range = (11.,13.5)
    group_m_range = (11.,16.)
    sat_frac_range = (0.,1.)
    Sigma_offset_range = (-0.01,0.01)
    
    # Number of bins in 2d plots (per axis)
    num_bins_2d = 20
    
    # Number of bins in 1d plots
    num_bins_1d = 50
    
    fig = plt.figure(4, figsize=(8,8))
    fig.subplots_adjust(hspace=0.001, wspace=0.001, left=0.08, bottom=0.095, top=0.975, right=0.93)
    
    # gridspec enables you to assign different formats to panels in one plot.
    gs = gridspec.GridSpec(4, 4, width_ratios=[1,1], height_ratios=[1,1])
    
    # Set up each of the 2D plots
    for  x_data,   x_range,        y_data,       y_range,     plot_shear, pos, bottom_row, left_col in (
        (sat_m,    sat_m_range,    group_m,      group_m_range,      True,  12, True,  True),
        (sat_m,    sat_m_range,    sat_frac,     sat_frac_range,     True,  13, True,  False),
        (sat_m,    sat_m_range,    Sigma_offset, Sigma_offset_range, False, 14, True,  False),
        (group_m,  group_m_range,  sat_frac,     sat_frac_range,     True,  8,  False, True),
        (group_m,  group_m_range,  Sigma_offset, Sigma_offset_range, False, 9,  False, True),
        (sat_frac, sat_frac_range, Sigma_offset, Sigma_offset_range, False, 4,  False, True)):
        
        # Work on the plot in this position
        plt.subplot(gs[pos])
        
        # Convert to 2d histogram for each of the shear, magnification, and overall data
        for index, label, colour in ((0, "Shear",         "b"),
                                     (1, "Magnification", "r"),
                                     (2, "Combined",      "k")):
            hist2D, xedges, yedges = np.histogram2d(x_data, y_data,
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
            cset = plt.contour(hist2D, levels, origin='lower',colors=[colour,colour],
                               linewidths=(2, 1),extent=extent,label=label)
            fmt = {}
            strs = [ '68\%', '95\%', ]
            for l,s in zip( cset.levels, strs ):
                fmt[l] = s
            plt.clabel(cset, inline=1, fontsize=10, fmt=fmt)
            for c in cset.collections:
                c.set_linestyle('solid')
        
        plt.xlim(x_range)
        plt.ylim(y_range)
        plt.ylabel(r'c',fontsize = 16)
        plt.xlabel(r'$\log_{\rm 10}{M_{\mathrm{sat,HDE}}}$ $(M_{\mathrm{\odot}})$',fontsize=16)
        ax.xaxis.set_major_locator(FixedLocator([11.3,11.4,11.5,11.6,11.7,11.8,11.9]))
        plt.yticks(fontsize=12)
        plt.xticks(fontsize=12)
        mass = np.linspace(11,12,100)
        conc = 4.67*(0.73*(10**mass)/(10**14))**(-0.11)
        _proxy = [plt.plot(0,0,'-',lw=3, color='red',linestyle='solid',label='HDE')]
        _proxy = [plt.plot(0,0,'-',lw=3, color='blue',linestyle='dashed',label='LDE')]
    plt.plot(mass,conc,'-',lw=3, color='black',linestyle='dotted',label='Standard M-c relation')
    plt.legend( loc='upper right',bbox_to_anchor=(-0.07, -0.3, 1, 1), bbox_transform=plt.gcf().transFigure)
    
    # Bin X,Y separately. As 1D bin, can use more bins now.
    S  = 26
    LSM = np.histogram(Satmass, bins=S, range=SMRANGE, normed=True)[0]
    LC = np.histogram(C, bins=S, range=CRANGE, normed=True)[0]
    # Restore positions lost by binning.
    SM = SMRANGE[0] + (SMRANGE[1]-SMRANGE[0])*np.array(range(0,len(LSM)))/float(len(LSM)-1)
    C = CRANGE[0] + (CRANGE[1]-CRANGE[0])*np.array(range(0,len(LC)))/float(len(LC)-1)
    LSM2 = np.histogram(Satmass2, bins=S, range=SMRANGE, normed=True)[0]
    LC2 = np.histogram(C2, bins=S, range=CRANGE, normed=True)[0]
    # Restore positions lost by binning.
    SM2 = SMRANGE[0] + (SMRANGE[1]-SMRANGE[0])*np.array(range(0,len(LSM2)))/float(len(LSM2)-1)
    C2 = CRANGE[0] + (CRANGE[1]-CRANGE[0])*np.array(range(0,len(LC2)))/float(len(LC2)-1)
    
    # bottom right panel: projected density of x.
    plt.subplot(gs[0])
    plt.plot(SM, LSM, '-', lw=3, color='red',linestyle='solid')
    plt.plot(SM2, LSM2, '-', lw=3, color='blue',linestyle='dashed')
    
    plt.yticks([])
    plt.xlim(SMRANGE)
    plt.ylim(0.0, 1.1*np.max(LSM2))
    ax = plt.gca()
    ax.xaxis.set_major_locator(FixedLocator([11.4,11.5,11.6,11.7,11.8,11.9]))
    plt.xticks(fontsize=12)
    plt.ylabel(r'$\log_{\rm 10}{M_{\mathrm{sat,HDE}}}$ $(M_{\mathrm{\odot}})$',fontsize=16)
    
    # top left panel: projected density of y.
    plt.subplot(gs[3])
    plt.plot(C, LC, '-', lw=3, color='red',linestyle='solid')
    plt.plot(C2, LC2, '-', lw=3, color='blue',linestyle='dashed')
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(FixedLocator([4,8,12,16,20]))
    plt.yticks([])
    plt.xticks(fontsize=12)
    plt.xlabel(r'c', fontsize=16)
    plt.ylim(0.0, 1.1*np.max(LC2))
    plt.xlim(CRANGE)
    
    plt.savefig('plot_contoursmc.png', format='png')
    plt.show()

if __name__ == "__main__":
    main(sys.argv)
