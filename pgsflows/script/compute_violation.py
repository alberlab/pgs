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
from alab.analysis import structuresummary
from subprocess import call
import json
import os
__author__  = "Hanjun Shin"
__credits__ = ["Nan Hua"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"
  
if __name__=='__main__':
	parser = argparse.ArgumentParser(description="compute_violation.py")
	parser.add_argument('-s', '--struct_dir', type=str, required=True)  
	parser.add_argument('-t', '--theta', type=str, required=True)    	
	parser.add_argument('-n', '--nstruct', type=int, required=True) 
	parser.add_argument('-o', '--violationFile', type=str, required=True) 
	parser.add_argument('-g', '--genome', type=str, required=True)
	parser.add_argument('-r', '--resolution', type=int, required=True)
	
	args = parser.parse_args()
	
	##################################################
	#summary
	##################################################
	print args.theta
	s = structuresummary(target=args.struct_dir, usegrp=args.theta, nstruct=int(args.nstruct) )
	    
	##################################################
    #violations
	##################################################
	data = {}
	if os.path.isfile(args.violationFile) :
		with open(args.violationFile, 'r') as file:
			data = json.load(file)
	
	data.update( {args.theta : {"totalViolations.mean" : s.totalViolations.mean(), 
								"totalRestraints.mean()" : s.totalRestraints.mean(), 
								"violation" : (s.totalViolations.mean()/s.totalRestraints.mean()) }
				})
	
	if data.has_key("pLast") :
		data.update( {"pLastActDist" : "from.%s.to.%s.actDist" % (data["pLast"], args.theta) })
	else :
		data.update( {"pLastActDist" : None})

	data.update( {"pLast" : args.theta} )
	
	with open(args.violationFile, 'w') as file:
		json.dump(data, file)