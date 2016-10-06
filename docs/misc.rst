Miscellaneous
=============


Organization of the .hmat file
------------------------------

PGS saves a probability matrix in a binary hdf5 file with extension hmat.
In the following we describe the file format. The file can be operated using alab api (see :doc:`alabapi <alab.contactmatrix>`)
    - matrix (n * n numpy 2D matrix storing the hic contacts)
    - idx (n * 4 numpy struct array consist of [chrom, start, end, flag] fields)
    - genome (cPickle serialized string)
    - resolution (cPickle serialized string)
    - applyedMethods (cPickle serialized string, storing all process done after creating the matrix)
    
Organization of the .hms file
-----------------------------

The output of each modeling step is dumped in a file with .hms extension. One .hms file contains information of a structure that evolve from the highest to lowest probability levels in theta list. Usually, only the final model extracted from this file is used for further analysis.
In the following we describe the file format. The file can be operated using alab api (see :doc:`alabapi <alab.modelstructures>`)
    - genome (cPickle serialized string)
    - resolution (cPickle serialized string)
    - idx (n * 4 numpy struct array consist of [chrom, start, end, flag] fields)
    - <group name 1> 
    - <group name 2>
    - ...

Group names are usually threshold values, each group consists the following components:
    - xyz (n * 3 array, coordinates)
    - r (n * 1 array, radii)
    - log (cPickled serialized string, all processing logs)
    - pym (serialized pym file content)

Organization of the .hss file
-----------------------------

A structure summary file, .hss can be generated once all final structures are successfully generated. 
In the following we describe the file format. The file can be operated using alab api (see :doc:`alabapi <alab.structuresummary>`)
    - genome (cPickle serialized string)
    - resolution (cPickle serialized string)
    - idx (n * 4 numpy struct array consist of [chrom, start, end, flag] fields)
    - usegrp (cPickle serialized string of threshold)
    - nbead (cPickle serialized int, number of representative beads in haploid)
    - nstruct (cPickle serialized int, number of total structures)
    - radius (n * 1 array of bead radii)
    - coordinates (nstruct * n * 3 array, coordinates for all structures)
    - score (1 * nstruct arry of optimization score)
    - consecutiveViolations
    - contactViolations
    - intraRestraints
    - interRestrains 
    
Activation distance file
------------------------

We start the modeling from random configurations for each structure and assign contacts between TADs that occur in 100% of the population (see `Methods <http://dx.doi.org/10.1073/pnas.1512577113>`_). Since we increase the number of restraints gradually, we need to assign contacts to a set of structures in the population at every A/M cycles. To do so, we compute the cummulative histogram and determine the distance cutoff according to the probability of a TAD pair in contact. The files are saved under ``result/actDist`` (the first 4 columns are TAD1, TAD2, contact probability, distance_cutoff_to_activate_contact).



