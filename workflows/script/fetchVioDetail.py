#!/usr/bin/env python
##
# Copyright (C) 2016 University of Southern California and
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

import sys
import numpy as np
import multiprocessing
import re
import argparse

import alab.files
import alab.matrix
import alab.utils
import alab.plots

__author__  = "Nan Hua"
__credits__ = ["Hanjun Shin"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"

def getVios(copystart,copyend,queue, struct_dir, prob):
    for i in range(copystart,copyend):
        sts     = alab.files.modelstructures(struct_dir+'/copy'+str(i)+'.hms',[prob])
        logs    = sts[-1].log
        vios    = re.findall(': (\d+.\d+) k = (\d+.\d+) (\d+.\d+)',logs)
        detail  = []
        for v in vios:
            dist,k,score = float(v[0]),float(v[1]),float(v[2])
            detail.append(1+(2*score/k)**0.5/dist)
        
        
        queue.put(detail)
    #-
    return 0

def listener(queue, prob, output_file):
    detail = []
    i=0
    while True:
        rs = queue.get()
        if rs == 'kill':
            break
        detail += rs
        i+=1
        #if i % 1000 == 0:
        #    print "%.3f %%" % (i*100.0/nstruct)
    #-
    detail = np.array(detail)
    #np.savetxt('%s_viodetail.txt'%(prob),detail,fmt='%f')
    v = np.percentile(detail,99) #99% percentile
    cutoff = np.percentile(detail,99.9)
    vio = detail[detail<cutoff]
    alab.plots.histogram(output_file, vio, 300,format='pdf',xlab='Violation Ratio',ylab='Frequency',histtype='stepfilled',color='g',line=v)
    
    return 0


def plotVio(prob, nstruct, struct_dir, output_file):
	pid = 10
	record = []
	manager = multiprocessing.Manager()
	queue = manager.Queue()
	watcher = multiprocessing.Process(target = listener,args=(queue, prob, output_file))
	watcher.start()

	for k in range(pid):
		start = k*(nstruct/pid)
		end   = (k+1)*(nstruct/pid)
		process = multiprocessing.Process(target=getVios,args=(start,end,queue, struct_dir, prob))
		process.start()
		record.append(process)

	for process in record:
		process.join()
    
	queue.put('kill')
	watcher.join()

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="report_summary.py")
	parser.add_argument('-p', '--prob', type=str, required=True)    #last frequency
	parser.add_argument('-n', '--nstruct', type=int, required=True) #current freq
	parser.add_argument('-s', '--struct_dir', type=str, required=True)  #probility matrix file in contactmatrix format
	parser.add_argument('-o', '--output_file', type=str, required=True) #current freq
	
	args = parser.parse_args()
	
	plotVio(args.prob, args.nstruct, args.struct_dir, args.output_file)
