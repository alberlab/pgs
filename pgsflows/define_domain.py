###### define_domain.py ######
#!/usr/bin/env python
"""
Pipeline description documentation.


DefineDomainFlow Arguments:
	Required Arguments:
		-i, --input    <UNDECIDED>
		-o, --output   <UNDECIDED>
	Optional Arguments:
		--flag	  <UNDECIDED>
"""
import argparse
import sys
import os

from pyflow import WorkflowRunner
from workflows.utils.args import add_pyflow_args
from workflows.utils.args import default_pyflow_args
from workflows.utils.args import extend_pyflow_docstring

#from workflows.utils.sys_utils import ensure_file_exists
#from workflows.config import SiteConfig

__version__ = "0.0.1"

class DefineDomainFlow(WorkflowRunner):
	def __init__(self, input, output):
		self.input = input
		self.output = output
			  

# input validation
		#ensure_file_exists(self.bam)

	def workflow(self):
		# args = [
				# 'command',
				# '-i %s' % self.bam
				# ]

		# if self.flag:
			# args.append('--flag')

		# self.addTask('task_id', ' '.join(args))
		print "DefineDomainFlow"
		
if __name__ == "__main__":
	parser = argparse.ArgumentParser(usage=extend_pyflow_docstring(__doc__))
	parser.add_argument('-i', '--input', type=str, required=True)
	parser.add_argument('-o', '--output', type=str, required=True)
	parser.add_argument('--flag', default=False, action="store_true")
	add_pyflow_args(parser)
	args = parser.parse_args()

	define_domain_wf = DefineDomainFlow(
			args.input,
			args.output
			)

	sys.exit(define_domain_wf.run(**default_pyflow_args(args)))