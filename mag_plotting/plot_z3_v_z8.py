#!/usr/bin/env python

import matplotlib
import os
import sys

import icebrgpy.columns_reader as cr
import matplotlib.pyplot as pyplot
import numpy as np
import subprocess as sbp


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


sqdeg_per_sqrad = (180/np.pi)**2

def main(argv):
    
    # Magic values
    
    data_table_name = "/home/brg/git/CFHTLenS_cat/Data/field_stats.dat"
    
    paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/Magnification_Method/"
    
    reader = cr.columns_reader()
    
    reader.add("z35_dens","lens_dens_z_0.3",factor=1./sqdeg_per_sqrad)
    reader.add("z85_dens","lens_dens_z_0.8",factor=1./sqdeg_per_sqrad)

    cols = reader.read(data_table_name)
        
    # Do the plotting now
        
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    
    ax.set_xlabel(r"Galaxies per deg$^2$ in $0.3 \leq z < 0.4$",labelpad=5,fontsize=16)
    ax.set_ylabel(r"Galaxies per deg$^2$ in $0.8 \leq z < 0.9$",labelpad=5,fontsize=16)
    
    ax.scatter(cols[reader.index("z35_dens")],cols[reader.index("z85_dens")])
            
    ax.legend(loc='lower right')
        
        
    # Save the figure
    outfile_name = os.path.splitext(data_table_name)[0] + "_z3_v_z8.eps"
    pyplot.writefig(outfile_name, format="eps", bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location
    cmd = "cp " + outfile_name + " " + paper_location
    sbp.call(cmd,shell=True)
            
    # Show the figure
        
    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)