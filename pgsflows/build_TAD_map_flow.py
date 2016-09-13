###### build_domain_contact_map.py ######
#
#!/usr/bin/env python
##
# Copyright (C) 2016 University of Southern California and
#                          Hanjun Shin
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

"""
Pipeline description documentation.


BuildTADMapFlow Arguments:
	Required Arguments:
		-c, --sys_config    
		-m, --matrixFile   
		-d, --domainFile
		-o, --outputFile
		
	Optional Arguments:
		--flag	  <UNDECIDED>
"""
import argparse
import sys
import os
import json
from subprocess import Popen, PIPE

from alab.pyflow import WorkflowRunner
from alab.args import add_pyflow_args
from alab.args import default_pyflow_args
from alab.args import extend_pyflow_docstring

#from workflows.utils.sys_utils import ensure_file_exists
#from workflows.config import SiteConfig

__version__ = "0.0.1"

class BuildTADMapFlow(WorkflowRunner):
	def __init__(self, input_config):
		self.input_config = input_config
			  
# input validation
		#ensure_file_exists(self.bam)

	def workflow(self):
		if type(self.input_config) == str :
			if not (os.path.exists(self.input_config) and os.access(self.input_config, os.R_OK) ):
				raise Exception('Cannot Find or Access input_config file %s' % self.input_config)
			with open(self.input_config, 'r') as data_file: 
				self.input_config = json.load(data_file)
		
		if not self.input_config.has_key('source_dir') :
			raise Exception('%s : Input config error, it does not have source_dir' % os.path.name(__file__))		
					
		if not self.input_config.has_key('input') :
			raise Exception('%s : Input config error, it does not have input' % os.path.name(__file__))
		else :
			if not self.input_config['input'].has_key('TAD_file') :
				raise Exception('%s : Input config error, it does not have TAD_file' % os.path.name(__file__))
			#if not self.input_config['input'].has_key('contact_map_file') :
			#	raise Exception('%s : Input config error, it does not have contact_map_file' % os.path.name(__file__))
		
		if not self.input_config.has_key('output_dir') :
			raise Exception('%s : Input config error, it does not have output_dir' % os.path.name(__file__))
				
		if not self.input_config.has_key('modeling_parameters') :
			raise Exception('%s : Input config error, it does not have modeling_parameters' % os.path.name(__file__))
			
		#else :
		#	if not self.input_config['modeling_parameters'].has_key('probMat') :
		#		raise Exception('%s : Input config error, it does not have probMat' % os.path.name(__file__))
		
		if not self.input_config.has_key('system') :
			raise Exception('%s : Input config error, it does not have system ' % os.path.name(__file__))
		else :
			if not self.input_config['system'].has_key('default_core') :
				raise Exception('%s : Input config error, it does not have default_core' % os.path.name(__file__))
			if not self.input_config['system'].has_key('max_memMB') :
				raise Exception('%s : Input config error, it does not have max_memMB' % os.path.name(__file__))
				
		
		
		domainFile = self.input_config['input']['TAD_file']
		nCores = self.input_config['system']['default_core']
		memMb = self.input_config['system']['max_memMB']
		genome = self.input_config['input']['genome']
		resolution = int( self.input_config['input']['resolution'] )
		outputfile = "%s/probMat/probMat.hdf5.hmat" % self.input_config['output_dir']
		if Popen("mkdir -p %s/probMat" % self.input_config['output_dir'], stderr=PIPE, stdout=PIPE, shell=True).wait() :
			raise Exception("Cannot create output_dir : %s/probMat" % self.input_config['output_dir'])
		
		python_path = Popen("which python", shell=True, stdout=PIPE).stdout.read().rstrip('\n')		
				
		matrixFile = None
		if self.input_config['input'].has_key('raw_matrix_file') :
			matrixFile = self.input_config['input']['raw_matrix_file']
			source = '%s/buildTADMap.py' % self.input_config['source_dir']
			
			args = [
				python_path,
				source, 
				'--matrixfile', matrixFile, 
				'--domainfile', domainFile,
				'--outputfile', outputfile,
				'--genome', genome, 
				'--resolution % i' % resolution
				]
		
		 	task_label = "buildTADMap_flow"
			self.addTask(label=task_label, command=' '.join(args), nCores=nCores, memMb=memMb, retryMax=3, retryWait=2, retryWindow=0, retryMode="all")	
			
		elif self.input_config['input'].has_key('prob_matrix_file_txt') :
			matrixFile = self.input_config['input']['prob_matrix_file_txt']
			source = '%s/hdf5_converter.py' % self.input_config['source_dir']
			
			args = [
				python_path,
				source, 
				'--matrixfile', matrixFile, 
				'--domainfile', domainFile,
				'--outputfile', outputfile,
				'--genome', genome
				]
		
		 	task_label = "hdf5ConvertFlow"
			self.addTask(label=task_label, command=' '.join(args), nCores=nCores, memMb=memMb, retryMax=3, retryWait=2, retryWindow=0, retryMode="all")	
			
		else :
			raise Exception('%s : Input config error, it does not have matrix_file, such as raw_matrix_file or prob_matrix_file_txt' % os.path.name(__file__))
				
			
		
		# python_path = Popen("which python", shell=True, stdout=PIPE).stdout.read().rstrip('\n')		
		# args = [
			# python_path,
			# source, 
			# '--matrixfile', matrixFile, 
			# '--domainfile', domainFile,
			# '--output_dir', output_dir,
			# '--genome', genome, 
			# '--resolution % i' % resolution
			# ]
		
		# task_label = "buildTADMap_flow"
		# self.addTask(label=task_label, command=' '.join(args), nCores=nCores, memMb=memMb, retryMax=3, retryWait=2, retryWindow=0, retryMode="all")		
		
			
			
		# source = '%s/buildTADMap.py' % self.input_config['source_dir']
		
		# matrixFile = self.input_config['input']['contact_map_file'] 
		# #if self.input_config['input']['contact_map_file_txt'] :
		# #	matrixFile = self.input_config['input']['contact_map_file_txt']
		# #elif self.input_config['input']['contact_map_file_hdf5'] :
		# #	matrixFile = self.input_config['input']['contact_map_file_hdf5']
			
		# domainFile = self.input_config['input']['TAD_file']
		# genome = self.input_config['input']['genome']
		# resolution = int( self.input_config['input']['resolution'] )
		
		# #outputFile = self.input_config['modeling_parameters']['probMat']
		# output_dir = "%s/probMat" % self.input_config['output_dir']
		# if Popen("mkdir -p %s" % output_dir, stderr=PIPE, stdout=PIPE, shell=True).wait() :
			# raise Exception("Cannot create output_dir : %s" % output_dir)
				
		# nCores = self.input_config['system']['default_core']
		# memMb = self.input_config['system']['max_memMB']
		
		# python_path = Popen("which python", shell=True, stdout=PIPE).stdout.read().rstrip('\n')		
		# args = [
			# python_path,
			# source, 
			# '--matrixfile', matrixFile, 
			# '--domainfile', domainFile,
			# '--output_dir', output_dir,
			# '--genome', genome, 
			# '--resolution % i' % resolution
			# ]
		
		# task_label = "buildTADMap_flow"
		# self.addTask(label=task_label, command=' '.join(args), nCores=nCores, memMb=memMb, retryMax=3, retryWait=2, retryWindow=0, retryMode="all")		
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(usage=extend_pyflow_docstring(__doc__))
	parser.add_argument('-i', '--input_config', type=str, required=True)
	#parser.add_argument('-m', '--matrixFile', type=str, required=True)
	#parser.add_argument('-d', '--domainFile', type=str, required=True)
	#parser.add_argument('-o', '--outputFile', type=str, required=True)

	#parser.add_argument('--flag', default=False, action="store_true")
	add_pyflow_args(parser)
	args = parser.parse_args()

	buildTADMapFlow_wf = BuildTADMapFlow(
			input_config = args.input_config
			)

	sys.exit(buildTADMapFlow_wf.run(**default_pyflow_args(args)))