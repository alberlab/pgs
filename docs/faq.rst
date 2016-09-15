Frequently Asked Questions
==========================


.. topic:: FAQ


    We are happy to update this list with your questions, please send inquiry to ``nhua@usc.edu``.


1. What PGS is for?
    It is a user-friendly software package to compute 3D genome structures from contact frequency data (e.g. Hi-C matrix). Since it generates a lot of structures (population), it is better to run it on HPC clusters. Be ready to give a generous amount of disk space (a typical intermediate structure file can be ~4MB. But once a structure summary file is created (~1.2 GB for 10,000 structures), the intermediate files can be deleted).


#. How long I expect the PGS to complete 1,000 structures?
    It depends on the computing power you assign it to. A typical M-step for 2 x 2000 TADs can take around one hour. It will also depend on the theta list you set (correspond to A/M iteration cycles). The lower theta value will give more restraints to optimize, thus the longer is an A/M cycle. If you have 1,000 cpus running for PGS, and there will be 10 A/M cycles, you might get the final population in ~10 hours.

#. Some nodes of my computing clusters crashed, how should I proceed PGS?
    No worries, PGS can resubmit the fail jobs for you and continues without problems.

#. I accidentally killed the terminal where PGS is run, what should I do?
    No worries, just go to the working directory and execute the PGS again using the exact same command (``sh runPgs.sh``). PGS is capable of tracking the last interruption and restarting the workflows smoothly.


#. PGS was terminated because of errors before creating any results, what should I do?
    In most cases the failures come from the input files, e.g. the matrix file or TAD file. 
    Here are some points to check:

        - Make sure all formating rules are met. 
        - The Hi-C matrix should represent a complete genome (include gaps so it is continues) in uniformly-sized bins.
        - TADs representation must be continues (include gaps and centromeres and label them as "gap" and "CEN", respectively. Every chromosome has one "CEN"). Of course, the size of TADs can differ.
        - All numeric must be integers (Hi-C counts, TAD genomic loci, etc.) except probability values (floating number from 0 to 1). 


#. What is TAD and how to get it for PGS run?
    TAD (Topologically Associating Domain) is a continuous genomic region within wich interact relatively frequently, whereas interactions across a TAD boundary occur relatively infrequently. Depending on the genome, their size can vary from tens of kb to a few Mb. We think this is a good chromosomal unit for 3D models. There are many TAD calling algorithm out there you can use. Ours is pretty much simple and quick, it's called `TopDom <https://doi.org/10.1093/nar/gkv1505>`_ and can be `downloaded here <http://zhoulab.usc.edu/TopDom>`_.


#. How many structures do I need to generate?
    We think it will depend on your analysis, but in general the more the better. The radial position of chromosomes or TADs, for example, are relatively stable so that hundreds to a thousand of structures are good enough. To reproduce the Hi-C map that comes from million of cells, we need around 10,000 structures. Beyond that we may gain marginal increase of correlations but computational costly. If the object of study, e.g. higher order interactions, are rare events like ~1%, consider getting at least 10,000 structures.


#. What is ``probMat.hdf5.hmat`` under ``result/probMat/`` folder?
    It is a TAD-TAD probability matrix obtained from your raw matrix or converted from your TAD-TAD probability matrix text file. This matrix can be used to get a replica population: copy it to a new directory, run the PGS-helper and choose option 3.




