#!/usr/bin/env python

# This is a sample script to use alab functions to summarize structure information from modeling information
#

import alab
import sys
import os.path

usage = sys.argv[0] + " <*.hms structure file dir> <# of structures> <probablility level> <output>\n\
hms file dir: model output file directory e.g. result/structures\n\
number of structures : the number of structures generated in PGS\n\
probabilility level : model theta threshold\n\
output : output filename, *.hss"

try:
    hmsdir  = sys.argv[1]
    nstruct = sys.argv[2]
    problvl = sys.argv[3]
    outfile = sys.argv[4]
except:
    print usage
    sys.exit()

ss = alab.structuresummary(hmsdir,problvl,nstruct,pid=10)

ss.save(outfile)

