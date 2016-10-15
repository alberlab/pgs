# PGS: Population-based Genome Structure modeling tools package

[![Build Status](https://travis-ci.org/alberlab/pgs.svg?branch=master)](https://travis-ci.org/alberlab/pgs)
[![Documentation Status](https://readthedocs.org/projects/pgs/badge/?version=latest)](http://pgs.readthedocs.io/en/latest/?badge=latest)

---

PGS is a population-based 3D genome-modeling package implemented in Python. 
The software takes Hi-C matrix and chromosome regions segmentation or topological associated domains(TADs) information, 
which then generates an ensemble of structure population. The software also automatically generates analysis reports, 
such as structure quality based on scoring parameters, plots of radial positions and contact frequency maps from the structures. 
The whole codes are wrapped in Python, and users can simply execute it one time. 

> Note:
> PGS software will run on high performance computing environment (HPC), such as sun grid engine and TORQUE (pbs script), as well as local machine. But, it is highly recommended to run on HPC because of computational resource and running time.

Here is the overview of PGS pipeline:

![PGS Overview] (docs/images/pgs_overview.png)

### Installation


Requirements:

- Python 2.7
- Python packages ``numpy``, ``scipy``, ``pandas``, ``h5py``, ``matplotlib`` , ``seaborn``
- IMP ([Integrative Modeling Package](https://integrativemodeling.org/))

Conda package is recommended to install all the requirements. Either [Anaconda](<https://www.continuum.io/downloads>) or 
the minimal [Miniconda](http://conda.pydata.org/miniconda.html) are suitable for managing required packages including IMP. If you use Miniconda, then you can install as follows:

```bash
    $ conda install numpy scipy pandas h5py matplotlib seaborn
```
Install IMP using conda:

```bash
    $ conda config --add channels salilab
    $ conda install imp
```
All other dependencies for imp and python packages will be automatically installed.

Then install PGS workflow packages:

```bash
    $ python setup.py install
``` 
### PGS Helper GUI


PGS package includes Graphical User Interface (GUI) based helper program for user to run pgs easily. 
User can generate command script (i.e. runPgs.sh) and configuration file(i.e. input_config.json) through the PGS Helper.


### Run PGS Helper

To initialize PGS Helper:

```bash
    $ java -jar PGSHelper.jar
```

The following GUI will appear:

![PGS Helper] (docs/images/pgs_helper.png)
   
### RUN PGS

User can run pgs package through the following command.

```
    $ PROJECT_DIR> sh runPgs.sh
``` 

### References

Tjong *et al.* `Population-based 3D genome structure analysis reveals driving forces in spatial genome organizations <http://dx.doi.org/10.1073/pnas.1512577113>`_. *PNAS* **113**, E1663-E1672 (2016).

