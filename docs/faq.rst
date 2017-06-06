Frequently Asked Questions
==========================


.. topic:: FAQ

    We are happy to update this list with your questions, please write them on the `GitHub page <https://github.com/alberlab/pgs/issues>`_.

#. What is PGS for?
    It is a user-friendly software package to compute 3D genome structures from contact frequency data (e.g. Hi-C matrix). Since it generates a lot of structures (population), it is better to run it on HPC clusters. Be ready to give a generous amount of disk space (a typical intermediate structure file can be ~4MB. But once a structure summary file is created (~1.2 GB for 10,000 structures), the intermediate files can be deleted).


#. Should I edit the text under "Optional argument list" of PGS helper?
    Yes, you should. Replace ``"qname"`` with your queue on HPC, but please do not delete the quote marks there (you may also delete this option and its value if you usually do not need to specify it when you submit jobs). Replace ``hh:mm:ss`` with number of hours, minutes, and seconds you wish to limit the time for a job to run. You can also add additional options with similar syntax (place a pair of quotes for each new option and its value, separated by a comma, and keep the brackets as it is).

#. I have matrix output from other pipeline, how can I convert my matrix file to acceptable input for PGS?
    We are adding support for the commonly used matrix format. However it would never be complete. There are various ways to convert your matrix to fit PGS input.
        - The easiest way is to output it as dense matrix format described in :doc:`quickstart`. 
        - For HiC-Pro users, one can use `HiC-Pro utilities tool <http://nservant.github.io/HiC-Pro/UTILS.html#hicpro2juicebox-sh>`_ to convert to juicebox \*.hic format.
        - For hiclib/mirnylib users, check out `Cooler package <https://github.com/mirnylab/cooler/>`_ CLI interface. Here is the example to convert hiclib output to cooler format:

        ::
            cooler cload --hiclib --assembly <genome_assembly> bins.bed your_fragment_dataset.hdf5 output_cooler.cool
    
        -
            
#. How long should I expect the PGS to complete 1,000 structures?
    It depends on the computing power you assign it to. A typical M-step for 2 x 2320 TADs can take around 45 minutes on our HPC cluster (2.6 GHz speed). It will also depend on the theta list you set (correspond to A/M iteration cycles). The lower theta value will give more restraints to optimize, thus the longer is an A/M cycle. If you have 1,000 cpus running for PGS, and there will be 10 A/M cycles, you might get the final population in ~8 hours.


#. Some nodes of my computing clusters crashed and some of PGS jobs were terminated, what should I do?
    No worries, PGS can resubmit the fail jobs for you automatically and continues without problems. If PGS is still running, you do not need to do anything, just wait.

    .. warning::  Do not alter ``pyflow.data/`` during PGS run. It contains logs and workflow state information. Deleting this folder will cause PGS to run from the beggining of the workflow again.


#. I accidentally killed the terminal where PGS is running, how should I proceed PGS?
    No worries, just go to the working directory and execute the PGS again using the exact same command (``sh runPgs.sh``). PGS is capable of tracking the last interruption and restarting the workflows from there without hasle (as long as the last workflow state recorded in ``pyflow.data/`` remains valid). 


#. PGS was terminated because of errors before creating any results, what should I do?
    You can first check the log files created under ``pyflow.data/logs/`` and try to fix that problems. However, in most cases the failures come from the input files, e.g. the matrix file or TAD file. 
    Here are some points to check while fixing the error(s):

        - Make sure all formating rules are met. 
        - The Hi-C matrix should represent a complete genome (include gaps so it is continues) in uniformly-sized bins.
        - TADs representation must be continues (include gaps and centromeres and label them as "gap" and "CEN", respectively. Every chromosome has one "CEN"). Of course, the size of TADs can differ.
        - All numeric must be integers (Hi-C counts, TAD genomic loci, etc.) except probability values (floating number from 0 to 1). 


#. What is TAD and how to get it for PGS run?
    TAD (Topologically Associating Domain) is a continuous genomic region within which interact relatively frequently, whereas interactions across a TAD boundary occur relatively infrequently. Depending on the genome, their size can vary from tens of kb to a few Mb. We think this is a good chromosomal unit for 3D models. There are many TAD calling algorithm out there you can use. Ours is pretty much simple and quick, it's called `TopDom <https://doi.org/10.1093/nar/gkv1505>`_ and can be `downloaded here <http://zhoulab.usc.edu/TopDom>`_.


#. How many structures do I need to generate?
    We think it will depend on your analysis, but in general the more the better. The radial position of chromosomes or TADs, for example, are relatively stable so that hundreds to a thousand of structures are good enough. To reproduce the Hi-C map that comes from million of cells, we need around 10,000 structures. Beyond that we may gain marginal increase of correlations but computational costly. If the object of study, e.g. higher order interactions, are rare events like ~1%, consider getting at least 10,000 structures.


#. What is ``probMat.hdf5.hmat`` under ``result/probMat/`` folder?
    It is a TAD-TAD probability matrix obtained from your raw matrix or converted from your TAD-TAD probability matrix text file. This matrix can be used to get a replica population: copy it to a new directory, run the PGS-helper and choose option 3.


#. Are the messages on screen while PGS is running saved somewhere?
    Yes, the stdout log messages are accumulated in your working directory under ``pyflow.data/logs/pyflow_log.txt``.


#. What is the ``pyflow.data/logs/pyflow_tasks_stdout_log.txt`` for?
    It contains detail information of all specific running jobs, i.e. processing matrix and modeling report (timing and scoring for all structures at all A/M cycles). For instance, you can search for "copy0.hms" in the log file and see how it performs from initial to final stages.



#. Any reference for the PGS?
    - Hua *et al.* `PGS: a dynamic and automated population-based genome structure software <https://doi.org/10.1101/103358>`_. *BioRxiv* (2017).
    
    - Kalhor *et al.* `Genome architectures revealed by tethered chromosome conformation capture and population-based modeling <http://dx.doi.org/10.1038/nbt.2057>`_. *Nat. Biotechnol.* **30**, 90-98 (2012).
    - Tjong *et al.* `Population-based 3D genome structure analysis reveals driving forces in spatial genome organizations <http://dx.doi.org/10.1073/pnas.1512577113>`_. *PNAS* **113**, E1663-E1672 (2016).






