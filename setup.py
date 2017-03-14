from distutils.core import setup

setup(
	name = 'pgs', 
	version = '1.0.1', 
	author = 'Nan Hua', 
	author_email = 'nhua@usc.edu', 
	url = '', 
	description = '3D Modeling Pipeline Workflows',
	packages=['alab'],
	package_data={'alab' : ['genomes/*']}
	)
