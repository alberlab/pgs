#!/usr/bin/env python

#Implementation of human TAD model in IMP 2.5

import alab.modeling
import alab.utils
import time
import IMP
import IMP.core
import sys
import re
import os
import numpy as np
import argparse
import json

__author__  = "Nan Hua"
__credits__ = ["Nan Hua","Ke Gong","Harianto Tjong", "Hanjun Shin"]

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"

#===================main entry========================
#def main(probfile,lastfb,currentfb, output_file, useLastCoordinates=False,actdistfile=None)
def main(input_config):

	if not input_config.has_key('modeling_parameters') :
			raise Exception('%s : Input config error, it does not have modeling_parameters' % os.path.name(__file__))
	else :
		if not input_config['modeling_parameters'].has_key('probMat') :
			raise Exception('%s : Input config error, it does not have probMat' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('last_theta') :
			raise Exception('%s : Input config error, it does not have last_theta' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('theta') :
			raise Exception('%s : Input config error, it does not have theta' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('output_file') :
			raise Exception('%s : Input config error, it does not have output_file' % os.path.name(__file__))
		if not input_config['modeling_parameters'].has_key('last_actDist') :
			raise Exception('%s : Input config error, it does not have last_actDist' % os.path.name(__file__))	
	
	probfile = str( input_config['modeling_parameters']['probMat'] )
	lastfb = str( input_config['modeling_parameters']['last_theta'] )
	currentfb = str( input_config['modeling_parameters']['theta'] )
	
	nucleusRadius = 5000.0
	if input_config['modeling_parameters'].has_key('nucleus_radius') :
		nucleusRadius = float( input_config['modeling_parameters']['nucleus_radius'] )
	
	chromosomeOccupancy = 0.2
	if input_config['modeling_parameters'].has_key('chr_occupancy') :
		chromosomeOccupancy = float( input_config['modeling_parameters']['chr_occupancy'] )
	
	output_file = str( input_config['modeling_parameters']['output_file'] )
	actdistfile = str( input_config['modeling_parameters']['last_actDist'] )
	
	useLastCoordinates = False
	#if currentfb != 'p100' :
	#	useLastCoordinates = True
	if actdistfile != None and os.path.isfile(actdistfile) :
		useLastCoordinates = True
		
	newmodel = alab.modeling.tadmodel(probfile,nucleusRadius=nucleusRadius,chromosomeOccupancy=chromosomeOccupancy,contactRange=1)
	tstart = time.time()

	if useLastCoordinates:
		lastxyz,r = alab.modeling.readCoordinates(output_file,lastfb)
	else:
		lastxyz = None
	newmodel.set_coordinates(lastxyz)

	#------adding restraints
	newmodel.set_excludedVolume()			  #1
	newmodel.set_nucleusEnvelope(kspring=1)	#2
	newmodel.set_consecutiveBeads(lowprob=0.6) #3


	if not useLastCoordinates:
		newmodel.CondenseChromosome(rrange=0.5)
		fmaxRestraints = newmodel.set_fmaxRestraints(kspring=1) #4
		totalIntra = 0
		totalInter = 0
		for rs in fmaxRestraints: #4
			if type(rs) == IMP.core.PairRestraint:
				totalIntra += 1
			else:
				totalInter += 1
		newmodel.updateScoringFunction()
		newmodel.SimulatedAnnealing(50000,1000,nc=2,nstep=1000)
		newmodel.mdstep_withChromosomeTerritory(500,2000)
		score = newmodel.cgstep(1000)
		newmodel.logger.debug(u"#of Intra restraints: %d #of Inter restraints: %d score: %.8f" %(totalIntra,totalInter,score))
		print "#of Intra restraints: %d #of Inter restraints: %d score: %.8f" %(totalIntra,totalInter,score)
	else:
		actdist = np.genfromtxt(actdistfile,usecols=(0,1,3),dtype=(int,int,float))
		intraContactRestraints, interContactRestraints = newmodel.set_contactRestraints(actdist)
		newmodel.updateScoringFunction()
		newmodel.SimulatedAnnealing(8000000,50000,nc=5, nstep=1000) #try1
		score = newmodel.cgstep(2000)
		newmodel.logger.debug(u"score try1: %.8f" % (score))
		print "score try1: %.8f" % (score)
		newmodel.SimulatedAnnealing(500000,5000,nc=8, nstep=1000) #try2
		score = newmodel.cgstep(2000)
		newmodel.logger.debug(u"score try2: %.8f" % (score))
		print "score try2: %.8f" % (score)
		newmodel.SimulatedAnnealing(100000,1000,nc=5, nstep=1000) #try3
		score = newmodel.cgstep(1000)
		newmodel.logger.debug(u"#of Intra restraints: %d #of Inter restraints: %d score: %.8f" % (len(intraContactRestraints),len(interContactRestraints),score))
		print "#of Intra restraints: %d #of Inter restraints: %d score: %.8f" % (len(intraContactRestraints),len(interContactRestraints),score)
		
	newmodel.logger.debug(u"===========  Optimization Starts ==========")
	print "===========  Optimization Starts =========="
	minscore = 10*newmodel.nbead
	if not useLastCoordinates:
		newmodel.updateScoringFunction()
		temp, score = newmodel.SimulatedAnnealing_Scored(100000,10000, nc=2, lowscore=newmodel.nbead)
		if score > minscore:
			temp, score = newmodel.SimulatedAnnealing_Scored(100000,10000, nc=2, lowscore=minscore)
	else:
		newmodel.shrinkingOptimization(drange=0.2,shrinkScore=20,minscore=minscore,interScale=0.1)

	newmodel.SimulatedAnnealing(100000,5000,nc=4, nstep=1500)
	newmodel.SimulatedAnnealing(10000,500,nc=4, nstep=1500)
	newmodel.SimulatedAnnealing(5000,500,nc=5, nstep=2500)

	newmodel.logger.debug(u"----- Relaxing Structure  -----")
	print "----- Relaxing Structure  -----"
	temp, score = newmodel.SimulatedAnnealing_Scored(5000,700,nc=4, nstep=1500, lowscore=10)
	if not useLastCoordinates:
		newmodel.mdstep_withChromosomeTerritory(500,1500)
		score=newmodel.cgstep(500)
		newmodel.mdstep_withChromosomeTerritory(300,1000)
		score1=newmodel.cgstep(500)
		newmodel.mdstep_withChromosomeTerritory(298,500)
		score=newmodel.cgstep(5000)
	else:
		newmodel.mdstep(500,1500)
		score=newmodel.cgstep(500)
		newmodel.mdstep(300,1000)
		score=newmodel.cgstep(500)
		newmodel.mdstep(298,500)
		score=newmodel.cgstep(5000)
		
	if np.isnan(score):
		f = open(copyname+'.txt','w')
		f.write(copyname)
		f.close()

	newmodel.logger.debug(u"Final score %.5f" %(score))
	print "Final score %.5f" %(score)
	newmodel.logger.debug(u"Total time spent is %.1fs" % (alab.utils.timespend(tstart)))
	print "Total time spent is %.1fs" % (alab.utils.timespend(tstart))
	newmodel.logger.debug(u"================================")
	print "================================"
	newmodel.logger.warning(u"Evaluating Contacts")
	print "Evaluating Contacts"
	newmodel.logger.warning(u"Consecutive Beads:")
	print "Consecutive Beads:"
	total = newmodel.evaluateRestraints(newmodel.consecutiveBeadRestraints)
	print "%d violations in total"%(total)
	newmodel.logger.warning(u"%d violations in total"%(total))
	newmodel.logger.warning(u"Contact Restraints:")
	print "Contact Restraints:"
	if useLastCoordinates:
		total = 0
		total += newmodel.evaluateRestraints(newmodel.intraContactRestraints)
		total += newmodel.evaluateRestraints(newmodel.interContactRestraints)
		print "%d violations in total"%(total)
		newmodel.logger.warning(u"%d violations in total"%(total))
	else:
		newmodel.evaluateRestraints(newmodel.fmaxRestraints)
	newmodel.saveCoordinates(output_file, currentfb)
	
	print "Simulation Finished."
	return 0

#=======================/main entry/===================

if __name__ == "__main__":
	
	# parser.add_argument('-p', '--probfile', type=str, required=True)  #probility matrix file in contactmatrix format
	# parser.add_argument('-l', '--lastfb', type=str, required=True)	#last frequency
	# parser.add_argument('-c', '--currentfb', type=str, required=True) #current theta
	# parser.add_argument('-o', '--output_file', type=str, required=True) # dir for saving structures, e.g. xxx/structures/model*/
	# parser.add_argument('-a', '--actdistfile', type=str, required=False, default=None) # the activation distance file by Astep.py e.g. model0.from.p010.to.p005.actDist
	
	#parser = argparse.ArgumentParser(description="MStep.py")
	#parser.add_argument('-i', '--input_config', required=True)  #probility matrix file in contactmatrix format
	#args = parser.parse_args()
	#main(args.input_config)
	
	#data = json.loads(args.input_config)
	#main(data)
	
	print sys.argv[1]
	data = json.loads( sys.argv[1] )
	main(data)
	
	# if re.search('p100',args.currentfb):
		# main(args.probfile,args.lastfb,args.currentfb,args.output_file,False,args.actdistfile)
	# else:
		# main(args.probfile,args.lastfb,args.currentfb,args.output_file,True,args.actdistfile)
 
