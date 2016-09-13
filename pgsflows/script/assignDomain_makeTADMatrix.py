#!/usr/bin/env python
"""
This is the current standard pipeline for HiC data processing
Using topdom output to generate domain level 
"""
import alab.matrix
import alab.files
import numpy as np
import sys

#format:
matrixfile = "GSE63525_GM12878_100kb_krnorm_rmDiag.hdf5"
domainfile = 'gm12878.100k.Ke.removepeaks.2320.ind'
outputfile = "GSE63525_GM12878_100kb_krnorm_rmDiag_TadLevel.hdf5"


#load matrix in to memory along with idx info, genome version, and resolution of the matrix
m = alab.matrix.contactmatrix(matrixfile)

#removeDiagonals 
#m.removeDiagonal()

#GenomeWide smoothing
#default setting: w=3,s=3,p=3,z=5
m.smoothGenomeWideHighValue(w=5,s=1,p=1,z=5)
m.save('GSE63525_GM12878_100kb_rmDiag.smoothedGenomeWide.hdf5')
m.makeIntraMatrix('chr1').plot('chr1_smoothed.png')
#load domain list information, skip the header and use 
#column 1 (chromosome); 
#column 2 (start site); 
#column 3 (end site)
#column 0 (suppose to be value for bedgraph file, here we use bead id)
#column 6 (flag)
#topdom header is as follows: bead_ID CHR from.coord to.coord from.bin to.bin chrbin.from chrbin.to tag"
domain = alab.files.bedgraph(domainfile,usecols=(1,2,3,0,6),skip_header=1)

#assign the matrix with the domain info
#use pattern to filter the domain list flags, pattern="domain" will ignore all regions that don't contain "domain" flag 
m.assignDomain(domain,pattern="")


#generate domain level matrix, using top 10% of the domain contact
newmatrix = m.iterativeFmaxScaling(23.4)
newmatrix.plot('heatmap_GSE63525_GM12878_100kb_rmDiag_TadLevel_median.png',log=False)
newmatrix.save(outputfile)
