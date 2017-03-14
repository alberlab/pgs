#!/usr/bin/env python
##
# Copyright (C) 2016 University of Southern California and
#                          Nan Hua
# 
# Authors: Nan Hua, Ke Gong, Harianto Tjong, and Hanjun Shin
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import alab.matrix
import alab.files
import numpy as np
import sys
import os
import argparse
import subprocess

__author__  = "Nan Hua"
__credits__ = ["Nan Hua","Ke Gong","Harianto Tjong", "Hanjun Shin"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"

#matrixfile = "GSE63525_GM12878_combined_100kb_raw.hdf5"
#domainfile = 'gm12878.100k.Ke.removepeaks.2320.ind'
#outputfile = "GSE63525_GM12878_100kb_krnorm_rmDiag_TadLevel.hdf5"
    
def main(matrixfile, domainfile, outputfile, genome, resolution) : #fileFormat):
    #load matrix in to memory along with idx info, genome version, and resolution of the matrix
    
#	if genome != 'hg19':
#		raise Exception('pgs v%s does not support %s genome yet' % (__version__, genome))
    
    #outputfile = '%s/probMat.hdf5.hmat' % output_dir
    if not os.path.isfile(matrixfile):
        raise IOError,"File %s doesn't exist!\n" % (matrixfile)
    if os.path.splitext(matrixfile)[1] == '.hic':
        m = alab.matrix.loadhic(matrixfile,genome=genome,resolution=resolution)
    elif os.path.splitext(matrixfile)[1] == '.cool':
        m = alab.matrix.loadcooler(matrixfile)
    else:
        m = alab.matrix.contactmatrix(matrixfile, genome=genome, resolution=resolution )
    #m = None
    #if fileFormat == 'hdf5' :
    #	m = alab.matrix.contactmatrix(matrixfile)
    #elif fileFormat == 'txt' :
    #	m = alab.matrix.contactmatrix(matrixfile, genome=genome, resolution=resolution )
    #else :
    #	raise Exception("File Format only accept either hdf5 or txt")
    
    #if you want to remove diagonal
    #use:
    m.removeDiagonal()

    #remove low coverage bins the lowest 1%
    m.removePoorRegions(cutoff=1)
    #bins that with high binomial split correlations are remained

    #identify Inter-chromosome contact outliers, run 100 iterations
    cutoff = m.identifyInterOutliersCutoff(N=100)
    #print cutoff

    #smooth inter-chromosome contact outliers using power law smoothing
    #default settings 
    #w=3,s=3,p=3
    m.smoothInterContactByCutoff(cutoff)

    #normalize the matrix using krnorm
    m.krnorm()
    #m.save(outputfile)
    #m.plot('heatmap_GSE63525_GM12878_100kb_krnorm.png')

    #format:
    #matrixfile = "GSE63525_GM12878_100kb_krnorm_rmDiag.hdf5"
    #domainfile = 'gm12878.100k.Ke.removepeaks.2320.ind'
    #outputfile = "GSE63525_GM12878_100kb_krnorm_rmDiag_TadLevel.hdf5"

    #load matrix in to memory along with idx info, genome version, and resolution of the matrix
    #m = alab.matrix.contactmatrix(matrixfile)

    #removeDiagonals 
    #m.removeDiagonal()

    #GenomeWide smoothing
    #default setting: w=3,s=3,p=3,z=5
    m.smoothGenomeWideHighValue(w=5,s=1,p=1,z=5)
    #m.save('GSE63525_GM12878_100kb_rmDiag.smoothedGenomeWide.hdf5') Hanjun
    #m.makeIntraMatrix('chr1').plot('chr1_smoothed.png') Hanjun
    
    #load domain list information, skip the header and use 
    #column 1 (chromosome); 
    #column 2 (start site); 
    #column 3 (end site)
    #column 0 (suppose to be value for bedgraph file, here we use bead id)
    #column 6 (flag)
    #topdom header is as follows: bead_ID CHR from.coord to.coord from.bin to.bin chrbin.from chrbin.to tag"
    #domain = alab.files.bedgraph(domainfile,usecols=(1,2,3,0,6),skip_header=1)
    # bed format
    domain = alab.files.bedgraph(domainfile,usecols=(0,1,2,2,3))

    #assign the matrix with the domain info
    #use pattern to filter the domain list flags, pattern="domain" will ignore all regions that don't contain "domain" flag 
    m.assignDomain(domain,pattern="")


    #generate domain level matrix, using top 10% of the domain contact
    newmatrix = m.iterativeFmaxScaling(25)
    #newmatrix.plot('heatmap_GSE63525_GM12878_100kb_rmDiag_TadLevel_median.png',log=False) Hanjun
    newmatrix.save(outputfile)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="buildTADMap.py")
    parser.add_argument('--matrixfile', type=str, required=True)  #raw matrix file
    parser.add_argument('--domainfile', type=str, required=True)  #domain file
    parser.add_argument('--outputfile', type=str, required=True)  #output dir
    parser.add_argument('--genome', type=str, required=True)      #genome	eg. hg19
    parser.add_argument('--resolution', type=int, required=True)  #resolution eg. 100000
    #parser.add_argument('--fileFormat', type=str, required=True)  #fileFormat hdf5|txt
        
    args = parser.parse_args()

    main(args.matrixfile, args.domainfile, args.outputfile, args.genome, args.resolution)
