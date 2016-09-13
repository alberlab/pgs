#!/usr/bin/env python

# Copyright (C) 2016 University of Southern California
# Frank Alber's lab
# Authors: Harianto Tjong & Nan Hua
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
__author__  = "HTjong & Nan Hua"
__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "htjong@usc.edu"


import numpy as np
import sys
import alab

def convertPDB(xyz,r,idx):
    """
    Given xyz r and index list, convert to pdb format string file content
    """
    assert len(idx)*2 == len(r)
    assert len(r) == len(xyz)
    #Locate centromere bead
    cenbead = np.flatnonzero(idx['flag']=='CEN')
    #output
    pdbtext = ""
    offset = -len(idx)
    for copy in ['A','B']:
        offset += len(idx)
        chrNum  = 0
        chrStart= 0
        for i in range(len(idx)):
            chrom = idx[i]['chrom']
            if chrom != idx[cenbead[chrNum]]['chrom']:
                chrNum += 1
                chrStart = i
            if i == cenbead[chrNum]:
                arm = "CEN"
            elif i < cenbead[chrNum]:
                arm = "PAM"
            else:
                arm = "QAM"
            
            chrchain = chr(chrNum+97) #convert to lowercase letters
            
            (x,y,z) = xyz[i+offset]
            line="ATOM %6d  %s %-3s %s %3d    %7.1f %7.1f %7.1f %5.0f\n"%(i+1,arm,copy+chrom.lstrip('chr'),chrchain,i-chrStart+1,x,y,z,r[i+offset])
            
            pdbtext += line
        #=
    #==
    return pdbtext

hmsfile = sys.argv[1]
problvl = sys.argv[2]
outfile = sys.argv[3]

hms = alab.modelstructures(hmsfile, [problvl])

pdb = convertPDB(hms[0].xyz,hms[0].r,hms[0].idx)
pdbfile = open(outfile,'w')
pdbfile.write(pdb)
pdbfile.flush()
pdbfile.close()



