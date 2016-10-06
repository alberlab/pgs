Tools
=====



.. topic:: Analysis tools to get results shown under the ``report/`` directory


        PGS package provides tools for analyzing the final structure population (see :doc:`alabapi`). The following are some examples of how to extract information from the structure-population.

        First, identify which model you would like to see, e.g. result/structure/copy0.hms. You can check group names in the file using a simple bash command:
	::

		$ h5ls result/structure/copy0.hms


	And something like this will appear:
            
	::

                0.01b                    Group
                0.01a                    Group
                0.01                     Group
                0.05                     Group
                0.1                      Group
                0.2                      Group
                1                        Group
                genome                   Dataset {SCALAR}
                idx                      Dataset {2320}


                                   

       	In that example, the "idx" has information of 2320 TADs and it saves iteration snapshots at theta = {1, .., 0.01, 0.01a, 0.01b}. Thus the final structure is in group "0.01b" (at theta level p=0.01 there are 3 A/M iteration cycles).


* Getting the 3D coordinates of the genome
	::

                import alab
                hmsfile = 'result/structure/copy0.hms
                problvl = '0.01a'
                hms = alab.modelstructures(hmsfile, [problvl])
                TADidx = hms.idx  #TADs information
                xyz = hms[0].xyz #diploid set of coordinates


    Now the user can use the coordinates, stored in ``xyz``, to do any analysis. The TAD information with genomic location is stored in ```TADidx``` variable. In the following we provide some other usage of coordinates.
	

* Getting the contact probability map 
    We show example lines on how to get contact probability maps under ``result/report/heatmap`` and ``result/report/intraMatrix```

	::

                import alab
                
                hmsfiledir = 'result/structure'
                problvl = '0.01a'
                nstruct = 1000
                summary = alab.structuresummary(hmsfiledir, problvl, nstruct)
                m = summary.getContactMap()
                m.plot('heatmap.png',format='png',clip_max=1)     
                m.makeIntraMatrix('chr1').plot('chr1_heatmap.pdf',format='pdf',clip_max=1)

* Getting the radial position
    Radial position of a TAD is calculated by the average radial positions across the structure population of that TAD. 0 marks the center of nucleus, and 1 marks nuclear envelope.

	::

                rp = summary.getBeadRadialPosition(beads=range(len(summary.idx)*2))
                rp_mean = rp.mean(axis=1)
                rp_hapmean = (rp_mean[:len(summary.idx)]+rp_mean[len(summary.idx):])/2

* Getting PDB
    Some users might wish to get the coordinates and radii in a PDB format, maybe for visualization purpose. Hence we provide some nice scripts under ``tool/`` directory. Simply execute the following shell command under ``$PROJECT_DIR/``:

	::

            $ tools/hms_export.py result/structure/copy0.hms 0.01b copy0.pdb

    The script takes 3 arguments (hmsfile, theta_group, and output_name), then a pdb file will be saved and it will look something like this:

	::

            ATOM      1  PAM A1  a   1     2899.1    58.6   855.0   218
            ATOM      2  PAM A1  a   2     3029.7   286.1  1257.0   244
            ATOM      3  PAM A1  a   3     2575.2   106.8  1117.7   202
            ....      .  ... ..  .   .     ......  ......  ......   ...
            ....      .  ... ..  .   .     ......  ......  ......   ...
            ATOM   1214  QAM BX  w  65    -2206.8   183.8  2465.6   202
            ATOM   1215  QAM BX  w  66    -2452.5   434.7  3049.5   238

..    Note::

    - The 2nd column marks the TADs ids.
    - PAM and QAM marks the short and long arms of a chromosome, respectively.
    - CEN marks the centromere representative TAD.
    - Chromosome homologue (4th column) is labeled as the chromosome name preceeded with A or B, e.g. A1 = the first homolog of chr1, BX = the 2nd homolog of chrX.
    - The first half of coordinates belong to the first diploid copy, the second half contains the homologues.
    - Chain name (5th column) is unique for each chromosome (e.g. chains "a" to "w" are for chr1 to chrX, respectively).
    - The 6th column contains TADs order in a chromosome.
    - Columns 6-8 record the 3D coordinates.
    - Column 9 stores the TADs radii.
    - In human genome 3D models, the nuclear radius is set to 5000 nm.



* Getting the structure-population summary file
    Instead of analyzing structures from many files, user can also aggregate the final structure population into a file and ignore the intermediate snapshots. Let the final population be in "001b" group, thus the following shell command should be run:

	::

            $ tools/hms_to_hss.py 0.01b structures_001b.hss

    This command outputs a summary file called "structures_001b.hss" which contains all coordinates of the last optimized structure population, their radii, TAD information, optimization scores, etc. At this point, if the user is not interested in the structures at intermediate steps, all ``structure/copy*.hms`` files can be deleted to release some disk space.

.. warning:: Check the content of the summary file (hss) first before deleting the \*.hms files!


