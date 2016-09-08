Quickstart
==========

Installation
------------

Requirements:

- Python 2.7
- Python packages ``numpy``, ``scipy``, ``pandas``, ``h5py``, ``matplotlib`` , ``six``
- IMP (`Integrative Modeling Package`_.)

.. _Integrative Modeling Package: https://integrativemodeling.org/

Conda package is recommended to install all the requirements. Either `Anaconda <https://www.continuum.io/downloads>`_ or 
the minimal `Miniconda <http://conda.pydata.org/miniconda.html>`_ are suitable for managing required packages including IMP. If you use Miniconda, then you can install as follows:

::

    $ conda install numpy scipy pandas h5py matplotlib six

Install IMP using conda:

::

    $ conda config --add channels salilab
    $ conda install imp

All other dependencies for imp and python packages will be automatically installed.

Then install PGS workflow packages:

::

    $ python setup.py install
    
PGS Helper GUI
--------------

PGS package includes Graphical User Interface (GUI) based helper program for user to run pgs easily. 
User can generate command script (i.e. runPgs.sh) and configuration file(i.e. input_config.json) through the PGS Helper.

.. Tip:: PGS Helper uses `Java Runtime Envrionment <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_, latest 8 update is recommended. 

Run PGS Helper
--------------

To initialize PGS Helper:

::

    $ java -jar PGSHelper.jar

The following GUI will appear:

.. image:: images/pgs_helper.png

