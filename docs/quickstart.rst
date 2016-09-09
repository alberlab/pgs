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
   :height: 1182px
   :width: 934px
   :scale: 50 %
   :align: center
   
A. PGS Project – Directory
~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify the project directory using the ``Browser`` button on right side. PGS will run in the specified project directory and all files, 
such as running script(i.e. ``runPGS.sh``), configuration file(i.e. ``input_config.json``),  log(i.e. ``pyflow.data``), and output results, 
will be stored in the directory.

B. PGS Source – Directory
~~~~~~~~~~~~~~~~~~~~~~~~~

Specify the directory of pgs source code using the ``Browser`` button on right side.

C. Input  
~~~~~~~~

1. Input Data

  *Option 1 : Raw(txt) + TAD*
  
     | Raw contact matrix file(txt) : First three columns are chromosome, start position(bp), and end position(bp) and followed by contact matrix. 
     | TAD_file(bed) : Chromatin Segmentation information or TAD file should be converted into `bed file format <https://genome.ucsc.edu/FAQ/FAQformat.html>`_. 
  
  *Option 2 : Prob(txt) + TAD*
  
     | Probability matrix file(txt)
     | TAD_file(bed)
     
  *Option 3 : Prob(hdf5)*
  
     | Probability matrix file(hdf5) : If the user have ever generated probability matrix using pgs, then user can use previous probability matrix. This process will skip the first workflow, buildTADMap task.

2. Genome : specify the genome version, Current PGS supports only hg19 with chromosome X.
3. Resolution : Specify the resolution of given input data 

D. Modeling Parameters
~~~~~~~~~~~~~~~~~~~~~~

1. Num of Structures : Specify the number of structures to generate using pgs. ``default = 1,000``
2. Violation Cutoff : Specify the violation cutoff.
   ``default = 0.05``
3. Theta Steps : Specify the list of thetas, 1 < theta < 0.
   ``default = 1, 0.2, 0.1, 0.05, 0.02, 0.01``
4. Max Iteration : Specify the number of maximum iterations for each theta.
   ``default = 10``

E. System Parameters
~~~~~~~~~~~~~~~~~~~~

1. Default Core : Specify the number of cores to use for default job, such as MStep.
2. Default MemMB : Specify the number of memory to use for default job, such as MStep. 
3. Max Core : Specify the number of cores to use for high demand jobs, such as AStep.
4. Max MemMB : Specify the number of memory to use for high demand jobs, such as AStep. 

F. Command Setup
~~~~~~~~~~~~~~~~

1. Run Mode : Specify the platform where pgs run on, such as Local, Sun Grid Engine or Torque. 
2. Core Limit : Specify the limit of number of cores to allow pgs to use based on user’s hpc policy.
3. Mem Limit : Specify the limit of memory to allow pgs to use based on user’s hpc policy.
4. Optional Argument List : Specify additional options for each job to run/be assigned correctly on user’s hpc, such as queue name and running time. Note that the option list will be applied to each job.
   For example: ``[‘-l’,’replace_your_qname_here’,’-l’,’walltime=333:00:00’]``

G. Generate Scripts 
~~~~~~~~~~~~~~~~~~~

Click ``Generate`` button on the bottom.



PGS Helper Output
-----------------

PGSInputGenerator creates input_config.json containing all input data address and parameters and running script (``runPGS.sh``) under the project directory. 

A. ``$PROJECT_DIR/input_config.json``

::

    {   "source_dir" : "[Directory name where pgs socurce is]",
        "input" : {
        "raw_matrix_file " : "[raw matrix file]",
            "TAD_file" : "[ TAD file, .bed format]"
            "resolution" : "[Resolution of input contact_map_file] e,g. 100000"
            "genome" : "[Genome version], e.g. hg19"
        },
        "output_dir" : "[Output Directory to store the results], e.g. $PROJECT_DIR/result",
        
        "modeling_parameters" : {
            "theta_list" : [Theta list] e.g, ["1", "0.2", "0.1","0.05","0.02","0.01"],
            "num_of_structures" : [Number of structure to generate] e.g. 1000,
            "max_iter_per_theta" : [Max Iterations per job] e.g. 10,
            "violation_cutoff" : [Violation Cutoff ] e.g. 0.05
        },
        "system" : {
            "max_core" : [Maximum number of cores in a single node], e.g. 8,
            "max_memMB" : [Maximum size of mem(MB) in a single node] e.g. 64000,
            "default_core" : [Default number of cores], e.g. 1,
            "default_memMB" : [Default size of mem(MB)] e.g. 1500
        }
    }

B. ``$PROJECT_DIR/runPGS.sh``

::

    python $PGS_DIRECTORY/pgs.py 
        --input_config $PROJECT_DIR/input_config.json 
        --run_mode [running platform] 
        --nCores 300 
        --memMb 800000 
        --pyflow_dir $PROJECT_DIR
        --schedulerArgList  ["-q","[qname]","-l","walltime=100:00:00"]

RUN PGS
-------

User can run pgs package through the following command.

::

    $ PROJECT_DIR> sh runPgs.sh
    
