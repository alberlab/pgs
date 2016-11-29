#!/usr/bin/env python

# Copyright (C) 2015 University of Southern California and
#                          Nan Hua
# 
# Authors: Nan Hua and Hanjun Shin
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

import argparse
import numpy as np
from subprocess import call

#from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import seaborn as sns
import alab.matrix
import sys
import os

sns.set(style="white")

from alab.analysis import structuresummary
from fetchVioDetail import plotVio
import alab.plots
import json

__author__  = "Nan Hua"
__credits__ = ["Hanjun Shin"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"
  
if __name__=='__main__':
	parser = argparse.ArgumentParser(description="report_summary.py")
	parser.add_argument('-s', '--struct_dir', type=str, required=True)  #probility matrix file in contactmatrix format
	#parser.add_argument('-p', '--prob', type=str, required=True)    #last frequency
	parser.add_argument('-n', '--nstruct', type=int, required=True) #current freq
	parser.add_argument('-o', '--output_dir', type=str, required=True) #current freq
	parser.add_argument('-g', '--genome', type=str, required=True)
	parser.add_argument('-r', '--resolution', type=int, required=True)
	parser.add_argument('-m', '--probMat', type=str, required=True)
	
	args = parser.parse_args()
	
	last_theta = 1
	violation_file = '%s/violation.json' % args.struct_dir
	if os.path.isfile( violation_file ) :
			with open(violation_file, 'r') as file:
				data = json.load(file)
	
			if data.has_key( "pLast" ) : 
				last_theta = data["pLast"]
			else :
				raise Exception("Cannot find violation rate for last_theta")
	else :
		raise Exception("Cannot find violation file, %s" % violation_file)
		
	
	chrList = []
	if args.genome == 'hg19':
		chrList = ['chr%i' % i for i in range(1,22)]
		chrList.append('chrX')
	elif args.genome == 'hg36':
		chrList = ['chr%i' % i for i in range(1,22)]
		chrList.append('chrX')
	elif args.genome == 'mm9' :
		chrList = ['chr%i' % i for i in range(1,19)]
		chrList.append('chrX')
	elif args.genome == 'mm10':
		chrList = ['chr%i' % i for i in range(1,19)]
		chrList.append('chrX')	
	else :
		raise Exception("Current %s genome is not supported" % args.genome)
	
	
	##################################################
	#summary
	##################################################
	#call(["mkdir", "-p", "%s/summary" % args.output_dir])
	s = structuresummary(target=args.struct_dir, usegrp=last_theta, nstruct=int(args.nstruct) )
	#s.save('%s/summary/summary.hss' % args.output_dir)
    
	##################################################
    #violations
	##################################################
	call(["mkdir", "-p", "%s/violations" % args.output_dir])
	with open("%s/violations/violations.txt" % args.output_dir, 'w') as file:
		file.write('Violation Ratio: %f' % (s.totalViolations.mean()/s.totalRestraints.mean()))
		#print 'Violation Ratio:',s.totalViolations.mean()/s.totalRestraints.mean()
		
	plotVio(last_theta, args.nstruct, args.struct_dir, "%s/violations/violation_plot.pdf" % args.output_dir)
	    
	##################################################
    #heatmap after modeling
	##################################################
	call(["mkdir", "-p", "%s/heatmap" % args.output_dir])
	m = s.getContactMap()
	m.plot('%s/heatmap/modeling_heatmap.pdf' % args.output_dir, format='pdf', clip_max=1)
	
	##################################################
	#IntraMatrix after modeling
	##################################################
	call(["mkdir", "-p", "%s/intraMatrix" % args.output_dir])
	for chr in chrList:
		m.makeIntraMatrix(chr).plot('%s/intraMatrix/%s_modeling_intraMatrix.pdf' % (args.output_dir, chr), format='pdf', clip_max=1)
		
	##################################################
    #heatmap of probMat
	##################################################
	call(["mkdir", "-p", "%s/heatmap" % args.output_dir])
	prob_mat = alab.matrix.contactmatrix(args.probMat)
	prob_mat.plot('%s/heatmap/probMat_heatmap.pdf' % args.output_dir, format='pdf', clip_max=1)
	
	##################################################
	#IntraMatrix of probMat
	##################################################
	call(["mkdir", "-p", "%s/intraMatrix" % args.output_dir])
	for chr in chrList:
		prob_mat.makeIntraMatrix(chr).plot('%s/intraMatrix/%s_probMat_intraMatrix.pdf' % (args.output_dir, chr), format='pdf', clip_max=1)
    	
	###################################################
    #radial position
	###################################################
	call(["mkdir", "-p", "%s/radialPlot" % args.output_dir])
	
	rp = s.getBeadRadialPosition(beads=range(len(s.idx)*2))
	rp_mean = rp.mean(axis=1)
	rp_hapmean = (rp_mean[:len(s.idx)]+rp_mean[len(s.idx):])/2
	
	np.savetxt('%s/radialPlot/radialPlot_summary.txt' % args.output_dir, rp_hapmean)
    
	for chr in chrList:
		s.plotRadialPosition('%s/radialPlot/%s_radialPlot.pdf' % (args.output_dir, chr), chr, format='pdf')
	
	
	###################################################
    #Compare Matrix
	###################################################
	#modelfile1 = "Ke_tadmodel30.hdf5"
	#modelfile2 = "Nan_tadmodel0.hdf5"
	cutoff = 0.05
	name1  = "Modeling"
	name2  = "Probability Matrix"
	comparePlotFile = "%s/heatmap/probMat_vs_modeling_compare_%i.pdf" % (args.output_dir, cutoff)

	#mx = alab.matrix.contactmatrix(modelfile1)
	mx = m
	my = prob_mat

	n = len(mx)
	x = mx.matrix[np.triu_indices(n,1)]
	y = my.matrix[np.triu_indices(n,1)]

	x1 = pd.Series(np.log(x[(x > cutoff) & (y > cutoff)]),name = name1)
	y1 = pd.Series(np.log(y[(x > cutoff) & (y > cutoff)]),name = name2)

	g = sns.jointplot(x1,y1,kind="kde",size=7,space=0)

	g.savefig(comparePlotFile)
            
        
    
    
