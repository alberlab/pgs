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

Specify the project directory using the ``Browse`` button on right side. PGS will run in the specified project directory and all files, 
such as running script(i.e. ``runPGS.sh``), configuration file(i.e. ``input_config.json``),  log(i.e. ``pyflow.data``), and output results, 
will be stored in the directory.

B. PGS Source – Directory
~~~~~~~~~~~~~~~~~~~~~~~~~

Specify the directory of pgs source code using the ``Browse`` button on right side.

C. Input  
~~~~~~~~

1. Input Data

  *Option 1 : Raw(txt) + TAD*
  
     | Raw contact matrix file (txt) : First three columns are chromosome, start position(bp), and end position(bp) and followed by contact matrix. 
     | TAD_file (bed) : Chromatin Segmentation information or TAD file should be converted into `bed file format <https://genome.ucsc.edu/FAQ/FAQformat.html>`. 
  
  *Option 2 : Prob(txt) + TAD*
  
     | Probability matrix file (txt)
     | TAD_file (bed)
     
  *Option 3 : Prob(hdf5)*
  
     | Probability matrix file (hdf5) : If the user have ever generated probability matrix using pgs, then user can use previous probability matrix. This process will skip the first workflow, buildTADMap task.

2. Genome : specify the genome version (current PGS supports only hg19 with chromosomes 1-22 and X).
3. Resolution : specify the resolution of given input data (in bp).

D. Modeling Parameters
~~~~~~~~~~~~~~~~~~~~~~

1. Num of structures : specify the number of structures to generate using pgs. ``default = 1,000``
2. Violation cutoff : specify the violation cutoff.
   ``default = 0.05``
3. Theta steps : specify the list of thetas, 1 < theta < 0.
   ``default = 1, 0.2, 0.1, 0.05, 0.02, 0.01``
4. Max iteration : specify the number of maximum iterations for each theta.
   ``default = 10``

E. System Parameters
~~~~~~~~~~~~~~~~~~~~
In order to proceed efficiently, PGS submits both single-core and multi-thread jobs on HPC clusters (e.g. for the M-step and A-step jobs, respectively).
Thus the following parameters need to be specified.
1. Default core : the number of cores to use for each regular job.
2. Default MemMB : the memory (Mb) to use for each regular job. 
3. Max cores : the number of cores to use for each multi-thread job.
4. Max MemMB : the total memory (Mb) to use for each multi-thread job. 

F. Command Setup
~~~~~~~~~~~~~~~~

1. Run mode : the platform where pgs run on, such as Local, Sun Grid Engine or Torque. 
2. Core limit : the maximum number of cores for PGS to use (limited to user’s quota).
3. Mem limit : the limit of memory for PGS to use.
4. Optional argument list : additional options for each job to run/be assigned properly on the user’s hpc, such as queue name, running time, etc. Note that the option list will be applied to each job.
   E.g. ``[‘-l’,’your_qname_here’,’-l’,’walltime=333:00:00’]``

G. Generate Scripts 
~~~~~~~~~~~~~~~~~~~

Click the ``Generate`` button on the bottom to write a file (input_config.json) with the parameters on the working directory which has been specified by the user.
There will be a confirmation window with ``Yes`` or ``No`` button, and at this point the user can see a simple instruction in the ``Usage`` box. If ``Yes`` is clicked, then the GUI will be closed.



PGS Helper Output
-----------------

PGSInputGenerator creates ``input_config.json`` containing all necessary information, and an execution script (``runPGS.sh``) under the project directory. 
At this point, the user tjust need to execute 


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
        --schedulerArgList  ["-q","qname","-l","walltime=100:00:00"]

RUN PGS
-------

User can run pgs package through the following command.

::

    $ PROJECT_DIR> sh runPgs.sh
    
