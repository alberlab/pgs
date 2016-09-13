from distutils.core import setup

setup(
	name = 'pgs', 
	version = '0.0.1', 
	author = 'Hanjun Shin', 
	author_email = 'shanjun@usc.edu', 
	url = '', 
	description = '3D Modeling Pipeline Workflows',
	packages=['alab'],
	package_data={'alab' : ['genomes/*']}
	)
