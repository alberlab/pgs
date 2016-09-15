Frequently Asked Questions
==========================


.. topic:: FAQ



1. How can I throw a question regarding PGS?
    We are happy to update this list with your questions, please send inquiry to ``nhua@usc.edu``.

#. How long I expect the PGS to complete 1,000 structures?
    It depends on the computing power you assign it to. A typical M-step for 2 x 2000 TADs can take around one hour. It will also depend on the theta list you set (correspond to A/M iteration cycles). The lower theta value will give more restraints to optimize, thus the longer is an A/M cycle. If you have 1,000 cpus running for PGS, and there will be 10 A/M cycles, you might get the final population in ~10 hours.

#. Some nodes of my computing clusters crashed, how should I proceed PGS?
    No worries, PGS can resubmit the fail jobs for you and continues without problems.

#. I accidentally killed the terminal where PGS is run, what should I do?
    No worries, just go to the working directory and execute the PGS again using the exact same command (``sh runPgs.sh``). PGS is capable of tracking the last interruption and restarting the workflows smoothly.


#. Can I rerun the PGS without using the raw matrix?
    Yes, you can. IF you want to run a replica calculation with the same input, just use the previously generated TAD-TAD probability matrix (``result/probMat/probMat.hdf5.hmat``), copy it to a new directory, run the PGS-helper and choose option 3.




