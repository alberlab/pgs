#!/usr/bin/env python

# Copyright (C) 2015 University of Southern California and
#                          Nan Hua
# 
# Authors: Nan Hua
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

__author__  = "Nan Hua"

__license__ = "GPL"
__version__ = "0.0.1"
__email__   = "nhua@usc.edu"

import os
import re
import math
import scipy
import time
import warnings
import numpy as np
from collections import namedtuple
from alabio import loadstream

#===========================================================================
class genome(object):
    def __init__(self,genomeName,usechr=['#','X']):
        datafile = os.path.join(os.path.dirname(os.path.abspath(__file__)) + '/genomes/' + genomeName + '.info')
        f = loadstream(datafile)
        self.info = np.genfromtxt(f,dtype=[('chrom','S5'),('length',int)])
        f.close() 
        removepattern = ''
        if not '#' in usechr: removepattern += '0-9'
        if not 'X' in usechr: removepattern += 'X'
        if not 'Y' in usechr: removepattern += 'Y'
        if not 'M' in usechr: removepattern += 'M'
        
        self.info = np.delete(self.info,np.nonzero([re.search('chr['+removepattern+']',c) for c in self.info['chrom']]))
    
    def bininfo(self,resolution):
        binSize    = [int(math.ceil(float(x)/resolution)) for x in self.info['length']]
        binStart   = [sum(binSize[:i]) for i in range(len(binSize))]
        
        chromList  = []
        binLabel   = []
        for i in range(len(self.info['chrom'])):
            chromList += [self.info['chrom'][i] for j in range(binSize[i])]
            binLabel     += [j for j in range(binSize[i])]
   
        startList  = [binLabel[j]*resolution for j in range(sum(binSize))]
        endList    = [binLabel[j]*resolution + resolution for j in range(sum(binSize))]
        
        binInfo    = [binSize,binStart,chromList,startList,endList]
        return namedtuple('binInfo', 'binSize, binStart, chromList, startList, endList')._make(binInfo)
  
    def getchrnum(self,chrom):
        findidx = np.flatnonzero(self.info['chrom']==chrom)
    
        if len(findidx) == 0:
            return -1
        else:
            return findidx[0]
  
    def getchrom(self,chromNum):
        assert isinstance(chromNum,int)
        return self.info['chrom'][chromNum]

#============================================end genome class
def timespend(t0):
    t1 = time.time()
    return int(round(t1-t0,0))

def listadd(a,b):
    """
        Add 2 array/list that have different sizes
    """
    if len(a) < len(b):
        c = np.array(b).copy()
        c[:len(a)] += a
    else:
        c = np.array(a).copy()
        c[:len(b)] += b
    return c

def boxplotStats(data):
    """
        The same as boxplot.stats() in R
        return 5 stats :
        Lower Fence
        Lower Quartile
        Median
        Upper Quartile
        Upper Fence
    """
    Q1,Q2,Q3 = np.percentile(data,[25,50,75])#get quartiles
    upperFence = Q3 + 1.5*(Q3-Q1)
    lowerFence = Q1 - 1.5*(Q3-Q1)
    return lowerFence,Q1,Q2,Q3,upperFence
  
def powerLawSmooth(matrix,target,w=3,s=3,p=3):
    """
        Power law smoothing function
        Given a matrix and a tuple (x,y), compute the smoothed value of (x,y)
        Parameters
        ----------
        matrix: numpy 2D array
        target: tuple of (x,y)
        w:      int of the window size, the smoothing is computed using target +/- w
        s:      weight of the location deviation
        p:      power of the location deviation
    """
    x,y = target
    csum = 0.0
    divider = 0.0
    for i in range(max(-w,-x),min(w+1,matrix.shape[0]-x)):
        for j in range(max(-w,-y),min(w+1,matrix.shape[1]-y)):
            decay = 1 / (abs(s*i) ** p + abs(s*j) ** p + 1.0)
            csum += matrix[x+i,y+j] * decay
            divider += decay
  
    return csum/divider

def smoothSpikesInBlock(matrix,w=3,s=3,p=3,z=5):
    """
        given a block matrix m, inspect whether any elements are above mean+z*sd
        replace them with the surrounding -w to +w local square window using power law smoothing
        Parameters:
        -----------
        matrix: block matrix in numpy.array
        w:      int of the window size, the smoothing is computed using target +/- w
        s:      weight of the location deviation
        p:      power of the location deviation
        z:      range of standard deviation
        
    """
    row,column     = matrix.shape
    print matrix.shape
    smoothedMatrix = np.copy(matrix)
    smoothedCounts = 0
    for i in range(row):
        for j in range(column):
            window  = matrix[max(i-w,0):min(i+w+1,row),max(j-w,0):min(j+w+1,column)].flatten()
            if matrix[i,j] > window.mean() + z*window.std(ddof=1):
                newVal = powerLawSmooth(matrix,(i,j),w,s,p)
                #newVal = powerLawSmooth(window,( i-max(i-w,0),j-max(j-w,0) ),w,s,p)
                if newVal < matrix[i,j]:
                    #print i,j,matrix[i,j],newVal
                    smoothedMatrix[i,j] = newVal
                    smoothedCounts += 1
            #--
        #--
    #--
    return smoothedMatrix,smoothedCounts

def binomialSplit(A,p=0.5):
    """
        split a matrix into 2 matrixes, using binomial randoms
    """
  
    assert isinstance(A,np.ndarray)
    split1 = np.zeros(np.shape(A))
    split2 = np.zeros(np.shape(A))
    if len(np.shape(A)) == 1:
        for i in range(np.shape(A)[0]):
            if A[i] > 0:
                s = np.random.binomial(A[i],p,1)[0]
                split1[i] = s
                split2[i] = A[i] - s
    elif len(np.shape(A)) == 2:
        for i in range(np.shape(A)[0]):
            for j in range(i,np.shape(A)[1]):
                if A[i,j] > 0:
                    s = np.random.binomial(A[i,j],p,1)[0]
                    split1[i,j] = s
                    split1[j,i] = s
                    split2[i,j] = A[i,j] - s
                    split2[j,i] = A[i,j] - s
    else:
        pass
    
    return split1,split2
     
def centerOfMass(xyz,r):
    """
    get center of mass of a list of particles, given xyz coordinates and radius
    """
    if len(xyz) != len(r):
        raise RuntimeError, "Dimension not agree"
    mass = r**3
    return np.sum(xyz*mass,axis=0)/sum(mass)

def intersectMatrixIndex(idx,chrom,querystart,querystop):
    """
    Fetch a intersection list with a given interval
    Return a list of indeces for with all records within the interval
    
    Parameter
    ---------
    idx : index structure array for alab.matrix.idx et al.
    chrom : chromsome e.g. 'chr1'
    querystart : start position
    queryend : end position
    
    Return
    ------
    array like indeces for all records 
    """
    import bisect
    
    intersectList = []
    sub = idx['chrom'] == chrom
    data = idx[sub]
    
    if len(data) > 0:
        rstart = np.flatnonzero(sub)[0]
        i = bisect.bisect(data['end'],querystart)
        while (i < len(data)) and (data[i]['start'] < querystop):
            intersectList.append(i)
            i += 1
    if len(intersectList) == 0:
        return None
    else:
        return np.array(intersectList) + rstart
    
def convertPDB(xyz,r,idx):
    """
    Given xyz r and index list, convert to pdb format string file content
    
    Parameter
    ---------
    xyz : numpy n*3 array, coordinates
    r : numpy n*1 array, radius
    idx : index structure array in alab.matrix.idx
    
    Return
    ------
    Converted pdb text string for future usage
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
