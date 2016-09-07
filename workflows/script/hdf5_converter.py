#!/usr/bin/env python

# Copyright (C) 2015 University of Southern California and
#                          Nan Hua
# 
# Authors: Hanjun Shin
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
from subprocess import call
import json
import os
__author__  = "Hanjun Shin"
__credits__ = ["Nan Hua"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"
  
if __name__=='__main__':
	parser = argparse.ArgumentParser(description="hdf5_converter.py")
	parser.add_argument('-i', '--matrixfile', type=str, required=True)  
	parser.add_argument('-d', '--domainfile', type=str, required=True)  
	parser.add_argument('-o', '--outputfile', type=str, required=True)    
	parser.add_argument('-g', '--genome', type=str, required=True) 
		
	args = parser.parse_args()
	import alab.matrix
	import alab.files
	# if matrixfile is txt file, write the output to output_file
	
	m      = alab.matrix.contactmatrix(args.matrixfile)
	domain = alab.files.bedgraph(args.domainfile,usecols=(0,1,2,2,3))
	m.assignDomain(domain)
	m._buildindex(m.domainIdx['chrom'],m.domainIdx['start'],m.domainIdx['end'],m.domainIdx['flag'])
	m.genome = args.genome
	m.resolution = 'TAD'
	m.save(args.outputfile)
	# if matrixfile is hdf5 file,
		# if output_file is not exist, then copy the matrixfile to output_file
		# if output_file is exist, nothing to do.
	
	