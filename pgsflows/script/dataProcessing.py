#!/usr/bin/env python
"""
This is the current standard pipeline for HiC data processing
Before topdom 
"""
import alab.matrix
import numpy as np
import sys

#format:
matrixfile = "GSE63525_GM12878_combined_100kb_raw.hdf5"
outputfile = "GSE63525_GM12878_100kb_krnorm_rmDiag.hdf5"

#load matrix in to memory along with idx info, genome version, and resolution of the matrix
m = alab.matrix.contactmatrix(matrixfile)

#if you want to remove diagonal
#use:
m.removeDiagonal()

#remove low coverage bins the lowest 1%
m.removePoorRegions(cutoff=1)
#bins that with high binomial split correlations are remained

#identify Inter-chromosome contact outliers, run 100 iterations
cutoff = m.identifyInterOutliersCutoff(N=100)
print cutoff

#smooth inter-chromosome contact outliers using power law smoothing
#default settings 
#w=3,s=3,p=3
m.smoothInterContactByCutoff(cutoff)

#normalize the matrix using krnorm
m.krnorm()
m.save(outputfile)
#m.plot('heatmap_GSE63525_GM12878_100kb_krnorm.png')