from distutils.core import setup

setup(
	name = 'workflows', 
	version = '0.0.1', 
	author = 'Hanjun Shin', 
	author_email = 'shanjun@usc.edu', 
	url = '', 
	description = '3D Modeling Pipeline Workflows',
	packages=['pyflow_alab', 'workflows', 'workflows.utils', 'workflows.modeling']
	)