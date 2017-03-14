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

GeneratePopulationOfGenomeStructure Arguments:
    Required Arguments:
        -i, --input_config    Path to intput configure file
"""
import argparse
import sys
import os
import os.path
import json
from subprocess import call

from alab.pyflow import WorkflowRunner
from alab.args import add_pyflow_args
from alab.args import default_pyflow_args
from alab.args import extend_pyflow_docstring

#from workflows.define_domain import DefineDomainFlow
from pgsflows.build_TAD_map_flow import BuildTADMapFlow
from pgsflows.modeling_structure_flow import ModelingStructureFlow
from pgsflows.report_flow import ReportFlow

#import workflows.build_TAD_map_flow.BuildTADMapFlow
#import workflows.modeling_structure_flow.ModelingStructureFlow
#import workflows.report_flow.ReportFlow
#from workflows.utils import ensure_file_exists

__version__ = "0.0.1"
__author__  = "Hanjun Shin"
__credits__ = ["Nan Hua"]

__license__ = "GPL"
__version__ = "1.0.1"
__email__   = "nhua@usc.edu"

class GeneratePopulationOfGenomeStructure(WorkflowRunner):
	def __init__(self, input_config) : #, contactMap, domainFile, output_dir, freq_list, nstruct):
		self.input_config = input_config
		
# input validation
        #ensure_file_exists(self.input_config)
		#ensure_file_exists(self.system_config)

	def workflow(self):
		if type(self.input_config) == str :
			if not (os.path.exists(self.input_config) and os.access(self.input_config, os.R_OK) ):
				raise Exception('Cannot Find or Access input_config file %s' % self.input_config)
			with open(self.input_config, 'r') as data_file: 
				self.input_config = json.load(data_file)
				
		if not self.input_config.has_key('source_dir') :
			raise Exception('%s : Input config error, it does not have source_dir key' % os.path.name(__file__))
		
		if not self.input_config.has_key('input') :
			raise Exception('%s : Input config error, it does not have input key' % os.path.name(__file__))
		else :
			if not self.input_config['input'].has_key('genome') :
				raise Exception('%s : Input config error, it does not have genome' % os.path.name(__file__))
			if not self.input_config['input'].has_key('resolution') :
				raise Exception('%s : Input config error, it does not have resolution' % os.path.name(__file__))
				
		if not self.input_config.has_key('modeling_parameters') :
			raise Exception('%s : Input config error, it does not have modeling_parameters key' % os.path.name(__file__))
			
		if not self.input_config.has_key('output_dir') :
			raise Exception('%s : Input config error, it does not have output_dir key' % os.path.name(__file__))
			
		if not self.input_config.has_key('system') :
			raise Exception('%s : Input config error, it does not have system key' % os.path.name(__file__))
										
		###############################################################
		# BuildTADMapFlow
		###############################################################
		build_task_config = self.input_config
		#build_task_config['source_dir'].update({'source_dir' : build_task_config['source_dir'] + "/workflows/script" })
		build_task_config['source_dir'] = build_task_config['source_dir'] + "/pgsflows/script" 
		
		#call(["mkdir", "-p", "%s/probMat" % build_task_config['output_dir']])
				
		build_TAD_map_wf = BuildTADMapFlow( build_task_config )
		build_TAD_map_task = "build_TAD_map_task"
		probMat = None
		
		if self.input_config['input'].has_key('raw_matrix_file') :
			if not self.input_config['input'].has_key('TAD_file') :
				raise Exception('%s : Input config error, it does not have TAD file' % os.path.name(__file__))
			else :
				self.addWorkflowTask(build_TAD_map_task, build_TAD_map_wf, dependencies=None)
				probMat = '%s/probMat/probMat.hdf5.hmat' % build_task_config['output_dir']
		elif self.input_config['input'].has_key('prob_matrix_file_txt') :
			if not self.input_config['input'].has_key('TAD_file') :
				raise Exception('%s : Input config error, it does not have TAD file' % os.path.name(__file__))
			else :
				self.addWorkflowTask(build_TAD_map_task, build_TAD_map_wf, dependencies=None)
				probMat = '%s/probMat/probMat.hdf5.hmat' % build_task_config['output_dir']
		elif self.input_config['input'].has_key('prob_matrix_file_hdf5') :
			build_TAD_map_task = None
			probMat = self.input_config['input']['prob_matrix_file_hdf5']
		else :
			raise Exception('%s : Input config error, it does not have matrix file' % os.path.name(__file__))
			
		
		####################################################
		#contact_map_file = None
		#if self.input_config['input'].has_key['contact_map_file'] :
		#	contact_map_file = self.input_config['input']['contact_map_file']
		#else :
		#	if not ( self.input_config['input'].has_key['prob_mat_file_txt'] or self.input_config['input'].has_key['prob_mat_file_hdf5'] ):
		#		raise Exception('%s : Input config error, you should give one of contact_map_file, prob_mat_file_txt, or prob_mat_file_hdf5' % os.path.name(__file__))
				
		#if contact_map_file == None :
		#	print "build_TAD_map_task is skipped, since it does not have contact_map_file(%s)" % contact_map_file
		#	build_TAD_map_task = None
		#else :
		#	self.addWorkflowTask(build_TAD_map_task, build_TAD_map_wf, dependencies=None)
		##################################################
		
		#prob_mat_file = None
		#if self.input_config['input'].has_key('prob_mat_file') :
		#	prob_mat_file = self.input_config['input']['prob_mat_file']
		#else :
		#	prob_mat_file = None
			
		##elif self.input_config['modeling_parameters'].has_key('probMat'). :
		##	prob_mat_file = self.input_config['modeling_parameters']['probMat']
		
		#if prob_mat_file == None :
		#	self.addWorkflowTask(build_TAD_map_task, build_TAD_map_wf, dependencies=None)
		#else :
		#	print "build_TAD_map_task is skipped, since it already has prob_mat_file(%s)" % prob_mat_file
		#	build_TAD_map_task = None
		
		###############################################################
		# ModelingStructureFlow
		###############################################################
		modeling_task_config = build_task_config
		#modeling_task_config['source_dir'].update({'source_dir' : modeling_task_config['source_dir'] + "/workflows/script" })
		#modeling_task_config['source_dir'] = modeling_task_config['source_dir'] + "/workflows/script" 
		
		#probMat = '%s/probMat/probMat.hdf5.hmat' % build_task_config['output_dir']
		#if prob_mat_file is not None :
		#	probMat = prob_mat_file
				
		modeling_task_config['modeling_parameters'].update({'probMat' : probMat})
		
		# modeling_structure_wf = ModelingStructureFlow(
			# sys_config = self.sys_config,
			# input = probMat,
			# output_dir = self.output_dir,
			# freq_list =  self.freq_list,
			# nstruct = self.nstruct
			# )
		modeling_structure_wf = ModelingStructureFlow(modeling_task_config)
		modeling_structure_task = "modeling_structure_task"
		self.addWorkflowTask(modeling_structure_task, modeling_structure_wf, dependencies=build_TAD_map_task)
		
		###############################################################
		# ReportFlow
		###############################################################
		report_task_config = modeling_task_config
		
		report_wf = ReportFlow(report_task_config)
		report_task = "report_task"
		self.addWorkflowTask(report_task, report_wf, dependencies=modeling_structure_task)
						
        #if self.flag:
        #    args.append('--flag')

        #self.addTask('task_id', ' '.join(args))
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(usage=extend_pyflow_docstring(__doc__))
	parser.add_argument('-i', '--input_config', type=str, required=True)
	#parser.add_argument('-c', '--sys_config', type=str, required=True)
	#parser.add_argument('-m', '--contactMap', type=str, required=True)
	#parser.add_argument('-d', '--domainFile', type=str, required=True)
	#parser.add_argument('-o', '--output_dir', type=str, required=True)
#    parser.add_argument('-s', '--system_config', default=True, action="store_true")
	add_pyflow_args(parser)
	args = parser.parse_args()

	GeneratePopulationOfGenomeStructure_wf = GeneratePopulationOfGenomeStructure(args.input_config)

	sys.exit(GeneratePopulationOfGenomeStructure_wf.run(**default_pyflow_args(args)))
