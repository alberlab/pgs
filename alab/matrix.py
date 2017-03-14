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
__credits__ = ["Nan Hua","Ke Gong","Harianto Tjong"]

__license__ = "GPL"
__version__ = "0.0.2"
__email__   = "nhua@usc.edu"

import numpy as np
import os.path
import re
import h5py
import copy
import cPickle
import warnings
from plots import plotxy, plotmatrix, histogram
import utils as alabutils

class contactmatrix(object):
    """
    A flexible matrix instant that supports various methods for processing HiC contacts
    
    Parameters
    ----------
    filename : str 
        matrix file stored in hdf5 format
        or an integer for the matrix size to initialize an empty matrix instance
    genome : str
        genome e.g.'hg19','mm9'
    resolution : int
        the resolution for the hic matrix e.g. 100000
    usechr : list
        containing the chromosomes used for generating the matrix
    
    Attributes
    ----------
    matrix : numpy 2d array 
        storing all infor for the hic contact matrix
    idx : numpy structure array 
        matrix index
    genome : str
        the genome
    resolution : int
        resolution for the contact matrix
    
    """
    _idxdtype = np.dtype([('chrom','S5'),('start',int),('end',int),('flag','S10')])
    def __init__(self,filename,genome=None,resolution=None,usechr=['#','X']):
        self._applyedMethods = {}
        if isinstance(filename,int):
            self.matrix=np.zeros((filename,filename),dtype = np.float32)
        elif isinstance(filename,str):
            if not os.path.isfile(filename):
                raise IOError,"File %s doesn't exist!\n" % (filename)
            if os.path.splitext(filename)[1] == '.hdf5' or os.path.splitext(filename)[1] == '.hmat':
                h5f = h5py.File(filename,'r')
                self.matrix = h5f['matrix'][:]
                self.idx    = h5f['idx'][:]
                if 'applyedMethods' in h5f.keys():
                    self._applyedMethods = cPickle.loads(h5f['applyedMethods'].value)
                
                if 'genome' in h5f.keys() and 'resolution' in h5f.keys():         
                    self.genome     = cPickle.loads(h5f['genome'].value)
                    self.resolution = cPickle.loads(h5f['resolution'].value)
                h5f.close()
            else:
                from alabio import loadstream
                f    = loadstream(filename)
                s    = f.next()
                line = re.split('\t+|\s+',s.rstrip())
                n    = len(line) - 3
                expectn = n
                if isinstance(genome,str) and isinstance(resolution,int):
                    genomedb    = alabutils.genome(genome,usechr=usechr)
                    bininfo     = genomedb.bininfo(resolution)
                    expectn     = len(bininfo.chromList)
                if expectn != n:
                    raise RuntimeError, "Dimension don't match, expected %s bins , get %s bins. Please check the input." %(expectn,n)
                idx  = []
                i    = 0
                tidx = line[0:3];tidx.append('')
                idx.append(tidx)
                self.matrix = np.zeros((n,n),dtype = np.float32)
                self.matrix[i] = line[3:]
                for s in f:
                    i += 1
                    line = re.split('\t+|\s+',s.rstrip())
                    tidx = line[0:3];tidx.append('')
                    idx.append(tidx)
                    self.matrix[i] = line[3:]
                f.close()
                self.idx    = np.core.records.fromarrays(np.array(idx).transpose(),dtype=self._idxdtype)
        else:
			raise RuntimeError, "Undefined input filename type!\n"
        #----------------end filename
        
        if isinstance(genome,str) and isinstance(resolution,int):
            if hasattr(self,"genome") and hasattr(self,"resolution"):
                raise RuntimeError, "Genome and resolution has already been specified."
            genomedb    = alabutils.genome(genome,usechr=usechr)
            bininfo     = genomedb.bininfo(resolution)
            flaglist    = ['' for i in range(len(bininfo.chromList))]
            self.genome = genome
            self.resolution = resolution
            self._buildindex(bininfo.chromList,bininfo.startList,bininfo.endList,flaglist)
      
    #==================================================
    def _buildindex(self,chromlist,startlist,endlist,flaglist):
        idxlist = np.column_stack([chromlist,startlist,endlist,flaglist])
        self.idx = np.core.records.fromarrays(np.array(idxlist).transpose(),dtype=self._idxdtype)
    def buildindex(self,**kwargs):
        warnings.warn("buildindex is deprecated, specify genome and resolution instead of building index manually.", DeprecationWarning)
        self._buildindex(**kwargs)
    #--------------------------------------------------
    def __str__(self):
        return self.matrix.__str__()
    def __repr__(self):
        return self.matrix.__repr__()
    def __len__(self):
        return self.matrix.__len__()

    def rowsum(self):
        return self.matrix.sum(axis=1)
    def columnsum(self):
        return self.matrix.sum(axis=0)
    def applyed(self,method):
        if method in self._applyedMethods:
            return True
        else:
            return False
    
    def _getZeroEntry(self):
        self.mask   = np.flatnonzero(self.rowsum() == 0)
    
    def _getMask(self,mask = None):
        if mask is None:
            self._getZeroEntry()
            return 0
        else:
            if isinstance(mask,np.ndarray):
                self.mask = mask
                return 1
            else:
                raise TypeError, "Invalid argument type, numpy.ndarray is required"
  
    def __checkGenomeResolution(self,genome,resolution):
        if hasattr(self,"genome") and hasattr(self,"resolution"):
            pass
        else:
            warnings.warn("No genome and resolution is specified within the file, try to assign attributes.")
            if (genome is None) or (resolution is None):
                raise ValueError, "No genome info is found! Genome and resolution parameter must be specified."
            else:
                self.genome = genome
                self.resolution = resolution

    #=========================================================filtering methods
    def removeDiagonal(self,force = False):
        if (not self.applyed('removeDiagonal')) or force:
            np.fill_diagonal(self.matrix,0)
            self._applyedMethods['removeDiagonal'] = True
        else:
            warnings.warn("Method removeDiagonal was done before, use force = True to overwrite it.")
      
    def removePoorRegions(self, cutoff=1, usepvalue = 0.1, force = False):
        """Removes "cutoff" percent of bins with least counts

        Parameters
        ----------
        cutoff : int, 0<cutoff<100
            Percent of lowest-counts bins to be removed
        usepvalue: float, 0<usepvalue<1
            use this pvalue as correlation cutoff to remove bins
            bins whose pvalue greater than this cutoff will be removed
        """
        if (not self.applyed('removePoorRegions')) or force:
            rowsum   = self.rowsum()
            mask     = np.flatnonzero(rowsum < np.percentile(rowsum[rowsum > 0],cutoff))
            #use correlation pvalue
            from scipy.stats import pearsonr,entropy
            newmask = []
            for i in mask:
                if rowsum[i] == 0:
                    continue
                split = alabutils.binomialSplit(self.matrix[i])
                corr,pvalue  = pearsonr(split[0],split[1])#returns corr
                rowentropy   = entropy(self.matrix[i])
                if pvalue > usepvalue:
                    newmask.append(i)
                    print i,rowsum[i],(corr,pvalue),"Remove",rowentropy
                else:
                    print i,rowsum[i],(corr,pvalue),"Keep",rowentropy
            self.matrix[newmask,:]    = 0
            self.matrix[:,newmask]    = 0
            self.idx['flag'][newmask] = 'Removed'
            print "%d low converage bins were removed." % (len(newmask))
            self._applyedMethods['removePoorRegions'] = (cutoff,len(newmask))
        else:
            warnings.warn("Method removePoorRegions was done before, use force = True to overwrite it.")
      
    def identifyInterOutliersCutoff(self,N=100):
        """
        Identify interchromosome outliers' cutoff
        Do an N round random choice as the original contact freq distribution and estimate the sample std for every contact freq
        If the sample std is larger than half of the frequency (contact #), lable this contact frequency as spourious
        cutoff is set to the number that first consecutive 2 non-spourious frequency from the right side (scan from high frequency to low)
        """
        if self.applyed('normalization'):
            raise RuntimeError, "Matrix is already normalized, raw matrix is needed."
        intermask = self.idx['chrom'][:,None] < self.idx['chrom'][None,:]
        interflatten = self.matrix[intermask]
        interflatten = interflatten[interflatten > 0]
        originHist = np.histogram(interflatten,interflatten.max())[0]
        repeatResults = np.zeros((N,interflatten.max()),dtype=int)
        for i in range(N):
            tmpChoice = np.random.choice(interflatten,len(interflatten))
            repeatResults[i] = alabutils.listadd(repeatResults[i],np.histogram(tmpChoice,tmpChoice.max())[0])
      
        comparison = np.std(repeatResults,axis=0) >= originHist/2
        i = len(comparison) -1
        while (comparison[i] or comparison[i-1]):i -= 1
        cutoff = i+1
        return cutoff
    
    def smoothInterContactByCutoff(self,cutoff,w=3,s=3,p=3,force=False):
        """
        given the cutoff, run a power law smoothing for the interchromosome matrix for contacts > cutoff
        """
        if (not self.applyed('smoothByCutoff')) or force:
            intermask = self.idx['chrom'][:,None] < self.idx['chrom'][None,:]
            x, y = np.nonzero(intermask)
            pos  = np.flatnonzero(self.matrix[intermask] > cutoff)
            for i in pos:
                cnew = alabutils.powerLawSmooth(self.matrix[ x[i]-w:x[i]+w+1 , y[i]-w:y[i]+w+1 ],(w,w),w,s,p)
                print x[i],y[i],self.matrix[x[i],y[i]],'-->',cnew
                self.matrix[x[i],y[i]] = cnew
                self.matrix[y[i],x[i]] = cnew
            self._applyedMethods['smoothByCutoff'] = cutoff
        else:
            warnings.warn("Method smoothInterContactByCutoff was done before,cutoff = %d use force = True to overwrite it."\
                            %(self._applyedMethods['smoothByCutoff']))
        
        
    #========================================================normalization methods
    def krnorm(self,mask = None,force=False,**kwargs):
        """
        using krnorm balacing the matrix (overwriting the matrix!)
        
        Parameters
        ----------
        
        mask : list/array 
            mask is a 1-D vector with the same length as the matrix where 1s specify the row/column to be ignored
            or a 1-D vector specifing the indexes of row/column to be ignored
            if no mask is given, row/column with rowsum==0 will be automatically detected and ignored
        large_mem : bool
            when large_mem is set to 1, matrix product is calculated using small chunks, 
            but this will slowdown the process a little bit.     
        """
        if (not self.applyed('normalization')) or force:
            from norm import bnewt
            self._getMask(mask)
            x = bnewt(self.matrix,mask=self.mask,check=0,**kwargs)*100
            self.matrix *= x 
            self.matrix *= x.T
            self._applyedMethods['normalization'] = 'krnorm'
        else:
            warnings.warn("Method %s was done before, use force = True to overwrite it." % (self._applyedMethods['normalization']))
  
    def vcnorm(self,iterations=1,mask = None,force=False):
        if (not self.applyed('normalization')) or force:
            self._getMask(mask)
            for i in range(iterations):
                print "\tIterations:",i+1
                rowsum   = self.rowsum()
                rowsum[self.mask] = 0
                totalsum = rowsum.sum()
                np.seterr(divide='ignore')
                rowsum   = 1/rowsum
                rowsum[self.mask] = 0
                self.matrix *= totalsum
                self.matrix *= rowsum 
                self.matrix *= rowsum.T
            self._applyedMethods['normalization'] = 'vcnorm'
        else:
            warnings.warn("Method %s was done before, use force = True to overwrite it." % (self._applyedMethods['normalization']))
    
    #-----------------------------------------------------------------------------
    
    def scale(self, cellaverage = 1):
        """
        Scale matrix so that average of cells is the given value. 
        By default, the rowsum will be the number of rows/columns
        """
        rowsum = self.rowsum()
        totalsum = rowsum.sum()
        try:
            self.mask
        except AttributeError:
            self._getMask()
        self.matrix = self.matrix / totalsum * (cellaverage * (len(rowsum)-len(self.mask)) * (len(rowsum)-len(self.mask)))
  
    def range(self,chrom):
        """
        return the index range for a give chromsome
        """
        rangeList = np.flatnonzero(self.idx['chrom'] == chrom)
        if len(rangeList)==0:
            raise ValueError, "%s is not found in the index" %(chrom)
        else:
            return (rangeList[0],rangeList[-1]+1)
  
    def makeIntraMatrix(self,chrom):
        """
        substract a chromsome matrix given a chromsome name
        
        Parameters
        ----------
        chrom : str, chromosome name e.g 'chr1'
        """
        if self.applyed('subMatrix'):
            warnings.warn("This is already a submatrix!")
        try:
            rstart,rend = self.range(chrom)
        except ValueError:
            raise ValueError, "%s is not found in the index. Possibly you are not using the genome wide matrix" %(chrom)
        submatrix   = contactmatrix(rend - rstart)
        submatrix.matrix = self.matrix[rstart:rend,rstart:rend]
        submatrix.idx    = np.core.records.fromrecords(self.idx[rstart:rend],dtype=self._idxdtype)
    
        if hasattr(self,"genome") and hasattr(self,"resolution"):
            submatrix.genome     = self.genome
            submatrix.resolution = self.resolution
        else:
            warnings.warn("No genome and resolution is specified, attributes are recommended for matrix.")
      
        submatrix._applyedMethods = copy.deepcopy(self._applyedMethods)
        submatrix._applyedMethods['subMatrix'] = chrom
        return submatrix
    
    def getICP(self,index):
        """
        return inter-chromosomal proportion of a given bin index
        """
        chrom = self.idx[index]['chrom']
        cstart,cend = self.range(chrom)
        intrasum = sum(self.matrix[index][cstart:cend])
        totalsum = sum(self.matrix[index])
        return 1-intrasum*1.0/totalsum
    
    #----------------------------------------------------------genome wide smoothing stuff:
    def smoothGenomeWideHighValue(self,w=3,s=3,p=3,z=5,force=False):
        """
        Use power law smoothing function to smooth high spikes in chromosomes blocks
        
        Parameters
        ----------
        w : int 
            the window size, the smoothing is computed using target +/- w
        s : int 
            weight of the location deviation
        p : int
            power of the location deviation
        z : int 
            range of standard deviation to set cutoff
        """
        if self.applyed('subMatrix'):
            raise RuntimeError, "This is a submatrix, genome wide smoothing cannot be applyed."
        if (not self.applyed('smoothGenomeWide')) or force:
            smoothed = 0
            chrlist = np.unique(self.idx['chrom'])
            if sum(np.diagonal(self.matrix)) < 1:
                v0 = np.diagonal(self.matrix,1)
                v1 = np.append(v0[0],v0)
                v2 = np.append(v0,v0[-1])
                np.fill_diagonal(self.matrix,v1+v2)
                #replace 0 with large number so most likely the neighborhood of diagonal are not smoothed
            for row in range(len(chrlist)):
                rstart,rend = self.range(chrlist[row])
                for column in range(row,len(chrlist)):
                    cstart,cend           = self.range(chrlist[column])
                    print "Smoothing block (%s,%s)" % (chrlist[row],chrlist[column])
                    tmpMatrix,smoothedCounts = alabutils.smoothSpikesInBlock(self.matrix[rstart:rend,cstart:cend],w,s,p,z)
                    self.matrix[rstart:rend,cstart:cend] = tmpMatrix
                    self.matrix[cstart:cend,rstart:rend] = tmpMatrix.T
                    if row == column:
                        smoothed += smoothedCounts
                    else:
                        smoothed += 2*smoothedCounts
            self._applyedMethods['smoothGenomeWide'] = (smoothed,"w=%d,s=%d,p=%d,z=%d" % (w,s,p,z))
            print "Genomewide smoothing finished, %d contacts smoothed" % (smoothed)
        else:
            warnings.warn("Method smoothGenomeWideHighValue was done before, %s %d values smoothed. use force = True to overwrite it."\
                            %(self._applyedMethods['smoothGenomeWide'][1],self._applyedMethods['smoothGenomeWide'][0]))
    #==============================================================Probabiliy matrix methods
    def getDomainMatrix(self,domainChrom,domainStartPos,domainEndPos,rowmask,minSize=1,maxSize=None):
        """
        Return a submatrix defined by domainChrom, domainStartPos, domainEndPos
        
        Parameters
        ----------
        domainChrom : str
            domain chromosome e.g. 'chr1'
        domainStartPos : int 
            start position e.g. 0
        domainEndPos : int 
            end position e.g. 700000
        minSize : int, > 0
            min domain size
        maxSize : int, optional
            max domain size, in bins
            if the domain is larger than a given number of bins, this function will return None
        """
        chrStartBin,chrEndBin = self.range(domainChrom)
        domainStartBin  = chrStartBin + int( domainStartPos/ float(self.resolution) )
        domainEndBin    = chrStartBin + int( np.ceil( domainEndPos/ float(self.resolution) ) )
        newmatrix = None
        if (domainEndBin - domainStartBin) >= minSize:
            if maxSize is None:
                newmatrix = self.matrix[domainStartBin:domainEndBin, domainStartBin:domainEndBin]
            elif (domainEndBin - domainStartBin) <= maxSize:
                newmatrix = self.matrix[domainStartBin:domainEndBin, domainStartBin:domainEndBin]
            else:
                return None
        else:
            return None
        maskloc = np.intersect1d(range(domainStartBin,domainEndBin),rowmask)
        maskloc = maskloc - domainStartBin
        newmatrix = np.delete(np.delete(newmatrix,maskloc,axis=0),maskloc,axis=1)
        return newmatrix
  
    def getfmax(self,method = 'UF',minSize=1,maxSize=2000,removeZero=False,boxplotTrim=False,offdiag=1,target='median'):
        """
        calculate fmax based on different methods
        
        Parameters
        ----------
        
        method : str
            NM #neighbouring max
            UF #uniform fmax
        target : str
            'mean'/'median'
        """
        if self.applyed('subMatrix'):
            raise RuntimeError, "This is a submatrix, genome wide fmax cannot be applyed."
        if self.applyed('probabilityMatrix'):
            raise RuntimeError, "This is already a probability matrix!"

        fmax = None
    
        if method == 'NM':#method neighbour max
            fmax = np.zeros(len(self))
            for chrom in np.unique(self.idx['chrom']):
                cstart, cend = self.range(chrom) #cend is increased by one
                for i in range(cstart+1,cend-1):
                    fmax[i] = min(self.matrix[i,i+1],self.matrix[i,i-1])
                fmax[cstart] = self.matrix[cstart,cstart+1] #p telomere
                fmax[cend-1] = self.matrix[cend-1,cend-2] #q telomere
        elif method == 'UF':#method uniform fmax
            if not hasattr(self,"domainIdx"):
                raise RuntimeError, "Please use assignDomain(domain_bedgraph,pattern) to assign domain INFO"
            print "Using minSize = %d, eliminating domains smaller than %dkb." % (minSize,minSize*self.resolution/1000)
            print "Using maxSize = %d, eliminating domains larger than %dkb." % (maxSize,maxSize*self.resolution/1000)
      
            #Get all intra domain interactions (upper triangle)
            upperTriangle = []
            skipDomains = 0
            print "Including Off Diagonal %d" %(offdiag)
            
            rowmask = np.flatnonzero(self.rowsum() == 0) #removed bins
            for domainRec in self.domainIdx:
                domainMatrix = self.getDomainMatrix(domainRec['chrom'],domainRec['start'],domainRec['end'],rowmask,minSize,maxSize)#domain matrix, eliminating domains that are larger than 20 bins
                if domainMatrix is None:
                    skipDomains += 1
                    continue
                upperTriangle.extend(domainMatrix[np.triu_indices(len(domainMatrix),offdiag)])#get the upper triangle
            #--------scaning finished
            print "%d domains are scanned, %d domains are eliminated." % (len(self.domainIdx),skipDomains)
            upperTriangle = np.array(upperTriangle)
            if removeZero:
                print "Removing zeros"
                upperTriangle = upperTriangle[upperTriangle > 0]
            lowerFence,Q1,Q2,Q3,upperFence = alabutils.boxplotStats(upperTriangle)#get quartiles and fence
            print lowerFence,Q1,Q2,Q3,upperFence
            if boxplotTrim:
                print "Trimming outliers"
                upperTriangle = upperTriangle[(upperTriangle > lowerFence) & (upperTriangle < upperFence)] #trim to better range
                lowerFence,Q1,Q2,Q3,upperFence = alabutils.boxplotStats(upperTriangle)#get quartiles and fence
                print lowerFence,Q1,Q2,Q3,upperFence
            if target == "median":
                fmax = Q2#get median
            elif target == "mean":
                fmax = upperTriangle.mean()#get mean
            else:
                raise RuntimeError, "target take only 'median' or 'mean' method"
        else:
            raise RuntimeError, "Please use legal method parameters:'NM' or 'UF'!"
        return fmax
    
    #-----------------------use fmax to get prob matrix
    def fmaximization(self,**kwargs):
        warnings.warn("fmaximization is deprecated, function name changed to fmaxScaling.", DeprecationWarning)
        self.fmaxScaling(**kwargs)
        
    def fmaxScaling(self,fmax,force=False):
        """
        use fmax to generate probability matrix
        for uniform fmax, simply divide the matrix by fmax and clip to 1
        for neighbouring contact fmax
        P[i,j] = F[i,j]/min(fmax[i],fmax[j])
        """
        if self.applyed('probabilityMatrix') and (not force):
            raise RuntimeError, "This is already a probability matrix!,use force to overwrite"
        if isinstance(fmax,float) or isinstance(fmax,np.float32) or isinstance(fmax,int):
            print "Uniform fmax detected"
            self.matrix = self.matrix/fmax
            self.matrix = self.matrix.clip(max=1)
            self._applyedMethods['probabilityMatrix'] = 'Uniform Fmax=%f' % (fmax)
        else:
            raise AttributeError, "Not supported fmax type!"
    
    def assignDomain(self,domain,pattern=''):
        """
        Load Domain information
        
        Parameters
        ----------
        
        domain : alab.files.bedgraph instance
            bedgraph for domain definition
        pattern : str 
            a string use to filter the flags in the bedgraph
        """
        from files import bedgraph
        if not isinstance(domain,bedgraph):
            raise TypeError,"Bedgraph instance required, see alab.files.bedgraph for more details"
        self.domainIdx = domain.filter(pattern)
    
    def _generateMedianSummaryMatrix(self,summaryBinStart,summaryBinEnd):
        N = len(summaryBinStart)
        X = np.empty((N,N),np.float32)
        for i in range(N):
            print "Filling X[%d] from A[%d] to A[%d]" % (i,summaryBinStart[i],summaryBinEnd[i]-1)
            istart = int(summaryBinStart[i])
            iend   = int(summaryBinEnd[i])
            for j in range(i,N):
                #print "Filling X[%d] from A[:,%d] to A[:,%d]" % (i,summaryBinStart[j],summaryBinEnd[j]-1)
                jstart = int(summaryBinStart[j])
                jend   = int(summaryBinEnd[j])
                submatrix = self.matrix[istart:iend,jstart:jend]
                #making sure that empty bins are removed
                out = np.nanmedian(submatrix)
                if np.isnan(out):
                    out = 0
                X[i,j] = X[j,i] = out
        return X
        
    def _generateTopMeanSummaryMatrix(self,summaryBinStart,summaryBinEnd,top=10,removeOutlier=True):
        N = len(summaryBinStart)
        X = np.empty((N,N),np.float32)
        for i in range(N):
            print "Filling X[%d] from A[%d] to A[%d]" % (i,summaryBinStart[i],summaryBinEnd[i]-1)
            istart = int(summaryBinStart[i])
            iend   = int(summaryBinEnd[i])
            for j in range(i,N):
                jstart = int(summaryBinStart[j])
                jend   = int(summaryBinEnd[j])
                submatrix = self.matrix[istart:iend,jstart:jend].flatten()
                submatrix = submatrix[~np.isnan(submatrix)] #get all non-nan value
                if len(submatrix) < 1:
                    out = 0
                else:
                    if removeOutlier:
                        lowerFence,Q1,Q2,Q3,upperFence = alabutils.boxplotStats(submatrix)
                        submatrix = submatrix[(submatrix>=lowerFence) & (submatrix<=upperFence)]
                    else:
                        pass
                    bound = np.percentile(submatrix,100-top)
                    out = np.mean(submatrix[submatrix>=bound])
                X[i,j] = X[j,i] = out
        return X
    
    def makeDomainLevelMatrix(self,method='topmean',top=10,removeOutlier=True):
        """
        Use domain INFO to generate Domain level matrix
        
        Parameters
        ----------
        
        method : str
            "topmean" or "median"
        top : int 0<top<100
            the top percentage to calculate the mean, top=10 means top 10% of the subdomain matrix
        removeOutlier : bool
            option to remove outlier using 1.5IQR
        """
        if self.applyed('domainLevel'):
            raise RuntimeError, "This is already a domain level matrix!"
        if not hasattr(self,"domainIdx"):
            raise RuntimeError, "Please use assignDomain(domain_bedgraph,pattern) to assign domain INFO"
            
        domainLevelMatrix = contactmatrix(len(self.domainIdx))
        domainLevelMatrix._buildindex(self.domainIdx['chrom'],self.domainIdx['start'],self.domainIdx['end'],self.domainIdx['flag'])
        domainLevelMatrix.genome = self.genome
        domainLevelMatrix.resolution = self.resolution
        domainLevelMatrix._applyedMethods = copy.deepcopy(self._applyedMethods)
        
        summaryBinStart = np.zeros(len(self.domainIdx))
        summaryBinEnd   = np.zeros(len(self.domainIdx))
        for i in range(len(self.domainIdx)):
            chrStartBin,chrEndBin = self.range(self.domainIdx[i]['chrom'])
            summaryBinStart[i]    = chrStartBin + int(self.domainIdx[i]['start'] / float(self.resolution))
            summaryBinEnd[i]      = chrStartBin + int(np.ceil(self.domainIdx[i]['end'] / float(self.resolution)))
        
        self._getMask()#removed bins in self.mask
        self.matrix[self.mask,:] = np.nan
        self.matrix[:,self.mask] = np.nan
        if method == 'topmean':
            domainLevelMatrix.matrix = self._generateTopMeanSummaryMatrix(summaryBinStart,summaryBinEnd,top,removeOutlier)
            domainLevelMatrix._applyedMethods['domainLevel'] = "%s/top=%d%s" % (method,top,'%')
        elif method == 'median':
            domainLevelMatrix.matrix = self._generateMedianSummaryMatrix(summaryBinStart,summaryBinEnd)
            domainLevelMatrix._applyedMethods['domainLevel'] = "%s" % (method)
        
        return domainLevelMatrix
    
    def iterativeFmaxScaling(self,domainAverageContacts=23.2,tol=0.01):
        """
        Automatic fmax scaling to get domain level matrix and match the rowsum average domain level matrix to 
        domainAverageContacts
        """
        fmax = self.getfmax()
        domainMean = 0
        originMatrix = copy.deepcopy(self.matrix)
        while abs(domainMean - domainAverageContacts)/domainAverageContacts > tol:
            print "fmax=%f"%(fmax)
            self.matrix = copy.deepcopy(originMatrix)
            self.fmaxScaling(fmax,force=True)
            domainLevelMatrix = self.makeDomainLevelMatrix()
            domainMean = domainLevelMatrix.rowsum().mean()
            fmax = fmax/domainAverageContacts*domainMean
        return domainLevelMatrix
        
    #==============================================================plotting methods
    def plot(self,figurename,log=False,**kwargs):
        """
        plot the matrix heat map
            
        Parameters
        ----------
        figurename : str
        log : bool
            if True, plot the log scale of the matrix
            if False, plot the original matrix
        clip_max : float
        clip_min : float
            2 options that will clip the matrix to certain value
        cmap : matplotlib.cm instance
            color map of the matrix
        label : str
            label of the figure
        """
        if log:
            plotmatrix(figurename,np.log(self.matrix),**kwargs)
        else:
            plotmatrix(figurename,self.matrix,**kwargs)
        
    def plotZeroCount(self,figurename,**kwargs):
        zeroCount = []
        for i in range(len(self.matrix)):
            zeros = len(np.flatnonzero(self.matrix[i] == 0))
            if zeros != len(self.matrix):
                zeroCount.append(zeros)
        #--endfor
        histogram(figurename,
                zeroCount,
                int(len(self.matrix)/100),
                xlab = '# of Zeros', ylab = 'Frequency',
                **kwargs)
        
    
    def plotSum(self,figurename,outlier=False,line=None,**kwargs):
        """
        Print the rowsum frequency histogram
        
        Parameters
        ----------
        figurename : string
            Name of the plot
        outlier: bool
            option to select plotting the outlier line, only functioning if 'line' parameter is set to None
        line : float/array/list
            draw vertical lines at a list of positions 
        """
        rowsum = self.rowsum()
        if line is None:
            if outlier:
                line = (np.percentile(rowsum,75) - np.percentile(rowsum,25))*1.5 + np.percentile(rowsum,75)
        
        histogram(figurename,
                rowsum[rowsum > 0],
                int(len(self.matrix)/100),
                xlab = 'Row sums', ylab = 'Frequency',
                line = line,
                **kwargs)
        
    #==============================================================save method
    def save(self, filename):
        """
            Save the matrix along with information in hdf5 file
        """
        if (filename[-5:] != '.hmat'):
            filename += '.hmat'
        h5f = h5py.File(filename, 'w')
        h5f.create_dataset('matrix', data=self.matrix, compression = 'gzip', compression_opts=9)
        h5f.create_dataset('idx', data=self.idx, compression = 'gzip', compression_opts=9)
        h5f.create_dataset('applyedMethods', data=cPickle.dumps(self._applyedMethods))
        if hasattr(self,"genome") and hasattr(self,"resolution"):
            h5f.create_dataset('genome',data = cPickle.dumps(self.genome))
            h5f.create_dataset('resolution',data = cPickle.dumps(self.resolution))
        else:
            warnings.warn("No genome and resolution is specified, attributes are recommended for matrix.")
        
        h5f.close()
#----------------------------------End of Class contactmatrix-------------------------------------------------
#==============================================================================================================

def loadh5dict(filename):
    h5f    = h5py.File(filename,'r')
    genome           = cPickle.loads(h5f['genome'].value)
    resolution       = cPickle.loads(h5f['resolution'].value)
    #genomeIdxToLabel = cPickle.loads(h5f['genomeIdxToLabel'].value)
    binNumber        = cPickle.loads(h5f['binNumber'].value)
    newMatrix = contactmatrix(binNumber,genome,resolution)
    newMatrix.matrix[:] = h5f['heatmap'][:]
    return newMatrix

def loadhic(filename,genome='hg19',resolution=100000,usechr=['#','X'],verbose=False):
    from . import straw
    
    tgenome = alabutils.genome(genome)
    bininfo = tgenome.bininfo(resolution)

    m = contactmatrix(len(bininfo.chromList),genome=genome,resolution=resolution,usechr=usechr)
    for chr1 in tgenome.info['chrom']:
        i = tgenome.getchrnum(chr1)
        for chr2 in tgenome.info['chrom']:
            j = tgenome.getchrnum(chr2)
            if i > j:
                continue
            if verbose:
                print chr1,chr2
            
            result = straw.straw("NONE",filename,chr1[3:],chr2[3:],'BP',resolution)
            for t in range(len(result[0])):
                x = int(result[0][t]/resolution) + bininfo.binStart[i]
                y = int(result[1][t]/resolution) + bininfo.binStart[j]
                m.matrix[x,y] = result[2][t]
                m.matrix[y,x] = result[2][t]
            #-
        #--
    #--
    
    return m

def loadcooler(filename,usechr=['#','X'],verbose=False):
    import scipy.sparse
    
    h5         = h5py.File(filename,'r')
    genome     = str(h5.attrs['genome-assembly'])
    resolution = h5.attrs['bin-size']
    nbins      = h5.attrs['nbins']
    
    if verbose:
        print genome, resolution
    
    tgenome = utils.genome(genome)
    bininfo = tgenome.bininfo(resolution)
    sp = scipy.sparse.csr_matrix((h5['pixels']['count'],h5['pixels']['bin2_id'],h5['indexes']['bin1_offset']),shape=(nbins,nbins))
   
    m = contactmatrix(len(bininfo.chromList),genome=genome,resolution=resolution,usechr=usechr)
    m.matrix = sp.toarray()[:len(m.idx),:len(m.idx)].astype(np.float32)
    h5.close()
    
    t = m.matrix.diagonal().copy()
    m.matrix[np.diag_indices(len(m.idx))] = 0
    m.matrix += m.matrix.T
    m.matrix[np.diag_indices(len(m.idx))] = t
    
    return m
