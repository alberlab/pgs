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

import re
import sys
import numpy as np
import alab.matrix
import alab.modeling
import multiprocessing
import os
import argparse
import json

__author__  = "Nan Hua"
__credits__ = ["Nan Hua","Ke Gong","Harianto Tjong", "Hanjun Shin"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"

#===========define functions
def readinCoordinates(copystart,copyend,lastfb,coor_shared,structdir,nstruct,nbead):
	arr = np.frombuffer(coor_shared.get_obj())
	modelcoor = arr.reshape((nstruct,nbead*2,3)) #10000 * 4640 * 3
	for i in range(copystart,copyend):
		try:
			xyz,r = alab.modeling.readCoordinates(structdir+'/copy'+str(i)+'.hms',lastfb)
			modelcoor[i][:] = xyz
		except RuntimeError:
			print "Can't find result for copy %s , %s" %(str(i),lastfb)
#--
def existingPortion(v, rsum):
	return sum(v<=rsum)*1.0/len(v)
def cleanProbability(pij,pexist):
	if pexist < 1:
		pclean = (pij-pexist)/(1.0-pexist)
	else:
		pclean = pij
	return max(0,pclean)

def calcActdist(jobqueue,fid,coor_shared,queue,probmat,nstruct,nbead,lastprob, r):
	arr = np.frombuffer(coor_shared.get_obj())
	modelcoor = arr.reshape((nstruct,nbead*2,3)) #10000 * 4640 * 3
	for job in jobqueue:
		i = job[0]
		j = job[1]
		pwish = probmat.matrix[i,j]
		plast = lastprob[i,j]
		
		d1 = np.linalg.norm(modelcoor[:,i,:] - modelcoor[:,j,:],axis=1) -r[i] - r[j] #L2 norm
		d2 = np.linalg.norm(modelcoor[:,i+nbead,:] - modelcoor[:,j+nbead,:],axis=1)-r[i]-r[j]
		if probmat.idx[i]['chrom'] == probmat.idx[j]['chrom']:
			dists = np.concatenate((d1,d2))
		else:
			d3 = np.linalg.norm(modelcoor[:,i,:] - modelcoor[:,j+nbead,:],axis=1) -r[i] - r[j]
			d4 = np.linalg.norm(modelcoor[:,i+nbead,:] - modelcoor[:,j,:],axis=1) -r[i] - r[j]
			dtable = np.column_stack((d1,d2,d3,d4))
			dists = np.sort(dtable,axis=1)[:,0:2].flatten() #fetch top 2 smallest distance
		#-
		sortdist = np.sort(dists)
		pnow = existingPortion(sortdist,r[i]+r[j])
		
		t = cleanProbability(pnow,plast)
		p = cleanProbability(pwish,t)
		if p>0:
			o = min(2*nstruct-1, int(round(2*p*nstruct)))
			res = '%4d %4d %5.3f %7.1f %5.3f %5.3f\n'%(i,j,pwish,sortdist[o],p,pnow)
			queue.put(res)
	return 0

def listener(queue,actFile):
	fout=open(actFile,'w')
	while True:
		res = queue.get()
		if res == 'kill':
			break
		fout.write(str(res))
		fout.flush()
	fout.close()
#===========end

#===================main entry========================
#def main(probfile,structdir,actFile,lastfb,currentfb,nstruct,pids,plastfile):
def main(input_config):

	if not input_config.has_key('modeling_parameters') :
			raise Exception('%s : Input config error, it does not have modeling_parameters' % os.path.name(__file__))
	else :
		if not input_config['modeling_parameters'].has_key('probMat') :
			raise Exception('%s : Input config error, it does not have probMat' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('struct_dir') :
			raise Exception('%s : Input config error, it does not have struct_dir' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('actDist') :
			raise Exception('%s : Input config error, it does not have actDist' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('last_theta') :
			raise Exception('%s : Input config error, it does not have last_theta' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('theta') :
			raise Exception('%s : Input config error, it does not have theta' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('num_of_structures') :
			raise Exception('%s : Input config error, it does not have num_of_structures' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('last_actDist') :
			raise Exception('%s : Input config error, it does not have last_actDist' % os.path.name(__file__))	
		if not input_config['modeling_parameters'].has_key('pids') :
			raise Exception('%s : Input config error, it does not have pids' % os.path.name(__file__))	

	probfile = str( input_config['modeling_parameters']['probMat'] )
	structdir = str( input_config['modeling_parameters']['struct_dir'] )
	actFile = str( input_config['modeling_parameters']['actDist'] )
	lastfb = str( input_config['modeling_parameters']['last_theta'] )
	currentfb = str( input_config['modeling_parameters']['theta'] )
	nstruct = int( input_config['modeling_parameters']['num_of_structures'] )
	pids = int( input_config['modeling_parameters']['pids'] )
	plastfile = input_config['modeling_parameters']['last_actDist']
	
	#################################
	#		parse arguments		#
	#################################
	probmat	= alab.matrix.contactmatrix(probfile)
	nbead	  = len(probmat)
	#targetfreq = float(re.search('\d+',currentfb).group(0))/100.
	
	getnum = re.compile(r'[^\d.]+')
	#targetfreq = float(getnum.sub('',currentfb))/100.
	targetfreq = float(getnum.sub('',currentfb))
	
	
	#################################
	#read in structure conformation #
	#################################

	xyz,r = alab.modeling.readCoordinates(structdir+'/copy0.hms',lastfb)#tempread one file
	coor_shared = multiprocessing.Array('d',nstruct*nbead*2*3) #10000 * 4640 * 3
	pid = 20
	readpool = []
	for k in range(pid):
		start = k*(nstruct/pid)
		end   = (k+1)*(nstruct/pid)
		process = multiprocessing.Process(target=readinCoordinates,args=(start,end,lastfb,coor_shared,structdir,nstruct,nbead))
		process.start()
		readpool.append(process)

	for process in readpool:
		process.join()
		
	###############################
	#read in last probability file#
	###############################
	
	lastprob = np.zeros((nbead,nbead))
	if not (plastfile == None or plastfile == 'null'):
		plastfile = str(plastfile)
		lastinfo = np.genfromtxt(plastfile)
		for i in range(len(lastinfo)):
			lastprob[int(lastinfo[i,0]),int(lastinfo[i,1])] = lastinfo[i,4]
		#-
	#--

	###############################
	#calculate activation distance#
	###############################


	jobqueue = []
	for i in range(pids):
		jobqueue.append([])
	q = 0
	for i in range(nbead):
		for j in range(i+1,nbead):
			if probmat.matrix[i,j] >= targetfreq:
				jobqueue[q].append((i,j))
				q+= 1
				if q >= pids:
					q = 0
				#-
			# skip those low frequency entries	

		#--
	#--
	record = []
	manager = multiprocessing.Manager()
	q = manager.Queue()
	watcher = multiprocessing.Process(target=listener,args=(q,actFile))
	watcher.start()

	for i in range(pids):
		process = multiprocessing.Process(target = calcActdist,args=(jobqueue[i],i,coor_shared,q,probmat,nstruct,nbead, lastprob,r))
		process.start()
		record.append(process)
	for process in record:
		process.join()
		
	q.put('kill')
	watcher.join()

	return 0
#=======================/main entry/===================

if __name__ == "__main__":
	
	# parser.add_argument('--probfile', type=str, required=True)  #probility matrix file in contactmatrix format <input>
	# parser.add_argument('--structdir', type=str, required=True) #xxx/structures/model*/ <input>
	# parser.add_argument('--actFile', type=str, required=True)   #activation distance file, <output>
	# parser.add_argument('--lastfb', type=str, required=True)	#last theta <parameter>
	# parser.add_argument('--currentfb', type=str, required=True) #current theta <parameter>
	# parser.add_argument('--nstruct', type=int, required=True)	#number of structures <parameter>
	
	# parser.add_argument('--pids', type=int, required=False, default=8)	#number of processers to use <setup>
	# parser.add_argument('--plastfile', type=str, required=False, default=None)#last activation distance file <optional input>
	#main(args.probfile, args.structdir, args.actFile, args.lastfb, args.currentfb, args.nstruct, int(args.pids), args.plastfile)
	
	#parser = argparse.ArgumentParser(description="AStep.py")
	#parser.add_argument('--input_config', type=str, required=True)
	#args = parser.parse_args()
	#main(args.input_config)
		
	data = json.loads( sys.argv[1] )
	main(data)
	