#!/usr/bin/env python

# This is a sample script to use alab functions to extract structure information from modeling information
#

import alab
import sys
import os.path

usage = sys.argv[0] + " <*.hms file> <probablility level> <output>\n\
hms file : model output file in result structures\n\
probabilility level : model theta threshold\n\
output : output filename, can be *.pym or *.pdb"

try:
    hmsfile = sys.argv[1]
    problvl = sys.argv[2]
    outfile = sys.argv[3]
except:
    print usage
    sys.exit()

hms = alab.modelstructures(hmsfile,[problvl])
if os.path.splitext(outfile)[1] == ".pym":
    hms[0].savepym(outfile)
elif os.path.splitext(outfile)[1] == ".pdb":
    hms[0].savepdb(outfile)
else:
    raise RuntimeError, "Output must be *.pym or *.pdb"
