###### Astep_flow.py ######
#!/usr/bin/env python
"""
Pipeline description documentation.

AStepFlow Arguments:

	Required Arguments:
		-i, --input_config 
		
"""
import argparse
import sys
import os
import json
from subprocess import Popen, PIPE
from subprocess import call

from alab.pyflow import WorkflowRunner
from alab.args import add_pyflow_args
from alab.args import default_pyflow_args
from alab.args import extend_pyflow_docstring

#from workflows.utils.sys_utils import ensure_file_exists
#from workflows.config import SiteConfig

__version__ = "0.0.1"
__AstepCores__ = 8

class AStepFlow(WorkflowRunner):
	#def __init__(self, sys_config, input, struct_dir, output, last_theta, theta, nstruct, last_actDist=None):
	def __init__(self, input_config):
		self.input_config = input_config
		# self.input = input
		# self.struct_dir = struct_dir
		# self.output = output
		# self.last_theta = last_theta
		# self.theta = theta
		# self.nstruct = nstruct
		
		# self.last_actDist = last_actDist		
		
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
					
		if not self.input_config.has_key('modeling_parameters') :
			raise Exception('%s : Input config error, it does not have modeling_parameters' % os.path.name(__file__))
		else :
			if not self.input_config['modeling_parameters'].has_key('probMat') :
				raise Exception('%s : Input config error, it does not have probMat' % os.path.name(__file__))
			if not self.input_config['modeling_parameters'].has_key('struct_dir') :
				raise Exception('%s : Input config error, it does not have structdir' % os.path.name(__file__))
			if not self.input_config['modeling_parameters'].has_key('actDist') :
				raise Exception('%s : Input config error, it does not have actDist' % os.path.name(__file__))
			if not self.input_config['modeling_parameters'].has_key('last_theta') :
				raise Exception('%s : Input config error, it does not have last_theta' % os.path.name(__file__))
			if not self.input_config['modeling_parameters'].has_key('theta') :
				raise Exception('%s : Input config error, it does not have theta' % os.path.name(__file__))	
			if not self.input_config['modeling_parameters'].has_key('num_of_structures') :
				raise Exception('%s : Input config error, it does not have num_of_structures' % os.path.name(__file__))
		
		if not self.input_config.has_key('system') :
			raise Exception('%s : Input config error, it does not have system ' % os.path.name(__file__))
		else :
			if not self.input_config['system'].has_key('max_core') :
				raise Exception('%s : Input config error, it does not have max_core' % os.path.name(__file__))
			if not self.input_config['system'].has_key('max_memMB') :
				raise Exception('%s : Input config error, it does not have max_memMB' % os.path.name(__file__))
		
		astep_src = '%s/Astep.py' % self.input_config['source_dir']
					
		#'/panfs/cmb-panasas2/shanjun/project/Hi-C/modelingPipe/workflows/AM/Astep.py',
		# args = [
				# 'python',
				# astep_src, 
				# '--probfile', self.input, 
				# '--structdir', self.struct_dir,
				# '--actFile', self.output,
				# '--lastfb', self.last_theta,
				# '--currentfb', self.theta,
				# '--nstruct', self.nstruct,
				# '--pids', '%d' % __AstepCores__
				# ]
				
		#if self.last_actDist != None :
		#	args.extend(['--plastfile', self.last_actDist])
		task_config = self.input_config
		task_config['modeling_parameters'].update({'pids' : self.input_config['system']['max_core']})
		
		# args = ['python',
				# astep_src,
				# '--input_config', json.dumps(task_config)]
				
		# args = ['python',
				# astep_src,
				# json.dumps(task_config)]
			
		#args = ['python',
		#		'/panfs/cmb-panasas2/shanjun/project/Hi-C/modelingPipe/workflows/test/test_code.py']

		python_path = Popen("which python", shell=True, stdout=PIPE).stdout.read().rstrip('\n')		
		task_label = "Astep_flow"
		self.addTask(label=task_label, command='%s %s \'%s\'' % (python_path, astep_src, json.dumps(task_config)), nCores=task_config['system']['max_core'], memMb=task_config['system']['max_memMB'], retryMax=3, retryWait=2, retryWindow=0, retryMode="all")
		#self.addTask(label=task_label, command=' '.join(args), nCores=1, memMb=1500, retryMax=3, retryWait=2, retryWindow=0, retryMode="all")
		# if self.flag:
			# args.append('--flag')
		
		#call(["mkdir", "-p", self.output_dir])
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(usage=extend_pyflow_docstring(__doc__))
	parser.add_argument('-i', '--input_config', type=str, required=True)
	
	# parser.add_argument('-i', '--input', type=str, required=True)
	# parser.add_argument('-s', '--struct_dir', type=str, required=True)
	# parser.add_argument('-o', '--output', type=str, required=True)
	# parser.add_argument('-l', '--last_theta', type=str, required=True)
	# parser.add_argument('-f', '--theta', type=str, required=True)
	# parser.add_argument('-n', '--nstruct', type=str, required=True)
		
	# parser.add_argument('-d', '--last_actDist', required=False, default=None)
	# parser.add_argument('--flag', default=False, action="store_true")
	add_pyflow_args(parser)
	args = parser.parse_args()

	# astep_wf = AStepFlow(
			# args.sys_config,
			# args.input,
			# args.struct_dir,
			# args.output,
			# args.last_theta,
			# args.theta,
			# args.nstruct,
			# args.last_actDist
			# )
	astep_wf = AStepFlow(args.input_config)

	sys.exit(astep_wf.run(**default_pyflow_args(args)))