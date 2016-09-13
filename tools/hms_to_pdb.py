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
import string
import sys
from itertools import izip
import re

import alab.matrix
import alab.analysis
import alab.utils


hmsfile = '../result/structure/copy0.hms'
problvl = ['0.1b']
outfile = 'copy0.pdb'

hms = alab.modelstructures(hmsfile, problvl)
bidx = hms.idx
nbead = bidx.shape[0]
chrs = list(set(bidx['chrom']))
nchr = len(chrs)
xyz = hms[0].xyz #diploid set of coordinates

#-----------------------------------------------------------------------
letters = string.letters[:26]
def pdbformat(n,a,ch,rid,resn,x,y,z,r):
    '''get a pdb formatted string '''
    line="ATOM %6d  %s %-3s %s %3d    %7.1f %7.1f %7.1f %5.0f\n"%(n,a,ch,rid,resn,x,y,z,r)
    return line

def chrblock(coor,rad,ids,pos,cen,chain,q):
    n = 0
    pdb_block = []
    for (x,y,z),r,i in izip(coor,rad,ids):
        n += 1
        if (i == cen):
            atmid = 'CEN'
        elif (i < q ):
            atmid = 'PAM'
        else:
            atmid = 'QAM'
        resid = pos + ch.lstrip('chr')
        pdb_block.append(pdbformat(i,atmid,resid,chain,n,x,y,z,r))
    return pdb_block

def com_rg(v,rd):
    mass = [l**3 for l in rd]
    xm1 = [(m1*i1,m1*j1,m1*k1) for (i1,j1,k1), m1 in izip(v,mass)]
    totm = sum(mass)
    comxyz = [sum(u)/totm for u in izip(*xm1)] #center of mass coordinates
    mma = sum( u**2 for u in comxyz )
    rcom = np.sqrt(mma)
    return rcom

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

chrs.sort(key=natural_keys)
#-----------------------------------------------------------------------

cenbead = np.flatnonzero(bidx['flag']=='CEN') #centromere representatives
cenbead_dict = {}
for ch,i in izip(bidx['chrom'][cenbead],cenbead):
    cenbead_dict[ch] = i


allcoord1 = xyz[:nbead] #first homologue
allcoord2 = xyz[nbead:] #second homologue
allrad = [hms[0].r[k][0] for k in range(nbead)]

fout=open(outfile,'w')
chrcount = 0 #for letters indicating chromosome
apdbs = []
bpdbs = []
for ch in chrs:
    chrchain = letters[chrcount]
    cen = cenbead_dict[ch]
    q = cen + 1
    chidx = np.flatnonzero(bidx['chrom']==ch)  #index or bead id, first column
    nb = len(chidx)
    chra_coor = [allcoord1[i] for i in chidx]
    chrb_coor = [allcoord2[i] for i in chidx]
    chr_rad = [allrad[i] for i in chidx]
    arcom = com_rg(chra_coor,chr_rad)
    brcom = com_rg(chrb_coor,chr_rad)
    if arcom < brcom:
        apdbs.append(chrblock(chra_coor,chr_rad,chidx,'A',cen,chrchain,q))
        bpdbs.append(chrblock(chrb_coor,chr_rad,chidx,'B',cen,chrchain,q))
    else:
        apdbs.append(chrblock(chrb_coor,chr_rad,chidx,'A',cen,chrchain,q))
        bpdbs.append(chrblock(chra_coor,chr_rad,chidx,'B',cen,chrchain,q))
    chrcount += 1
for line in apdbs:
    for l in line:
        fout.write(l)
for line in bpdbs:
    for l in line:
        fout.write(l)
fout.close()

sys.exit()
