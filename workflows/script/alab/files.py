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

import numpy as np
import re
import os
import bisect
import h5py
import cPickle

#====================================================================================
class bedgraph(object):
    """
    Required fields

        chrom      - name of the chromosome or scaffold. Any valid seq_region_name can be used, and chromosome names can be given with or without the 'chr' prefix.
        chromStart - Start position of the feature in standard chromosomal coordinates (i.e. first base is 0).
        chromEnd   - End position of the feature in standard chromosomal coordinates
        dataValue  - Track data values can be integer or real
    """
    _bedgraphdtype = np.dtype([('chrom','S5'),('start',int),('end',int),('value',float),('flag','S20')])
  
    def __init__(self,filename=None,usecols=(0,1,2,3,3),**kwargs):
        self.data            = {}
        self.__sorted_keys   = []
        self.itr             = 0
        if not filename is None:
            if isinstance(filename,str):
                from alabio import loadstream
                f = loadstream(filename)
                readtable = np.genfromtxt(
                                    f,
                                    dtype=self._bedgraphdtype,
                                    usecols=usecols,
                                    **kwargs)
                f.close()
            elif isinstance(filename,np.ndarray) or isinstance(filename,list):
                readtable = np.core.records.fromarrays(
                                    np.array(filename).transpose()[[usecols]],
                                    dtype = self._bedgraphdtype)
            for line in readtable:
                chrom = line['chrom']
                #if np.isnan(line['value']):
                    #line['value'] = 0
                if not chrom in self.data:
                    self.data[chrom] = []
                    self.data[chrom].append(line)
                else:
                    #self.data[chrom] = np.append(self.data[chrom],line)
                    self.data[chrom].append(line)

            for chrom in self.data:
                self.data[chrom] = np.core.records.fromrecords(self.data[chrom],
                                                               dtype = self._bedgraphdtype)
                self.data[chrom].sort(kind='heapsort',order='start')
      
            self._flush()
      
    #========================================================  
    def genchrnum(self,chrom):
        """ Sort by chromosome """
        if chrom:
            num = chrom[3:]
            if   num == 'X': num = 23
            elif num == 'Y': num = 24
            elif num == 'M': num = 25
            else: num = int(num)
        else:
            num = 0
        return num
  
    def _flush(self):
        self.__sorted_keys = sorted(self.data.keys(),key=lambda x:self.genchrnum(x))
        
    def __repr__(self):
        represent = ''
        for chrom in self.__sorted_keys:
            represent += '\n'+chrom+'\twith '+str(len(self.data[chrom]))+' records'
        return represent
  
    def __len__(self):
        length = 0
        for chrom in self.data:
            length += len(self.data[chrom])
        return length
  
    def __getonerec(self,key):
        """For cases a[i], output ith record"""
        if key < 0:
            key += len(self)
        if key > len(self):
            raise IndexError, "The index (%d) is out of range" % key
        for chrom in self.__sorted_keys:
            if key+1 - len(self.data[chrom]) > 0:
                key = key - len(self.data[chrom])
            else:
                return self.data[chrom][key]
      
    #++++++++++++++++++++++++++++++++++++++++++++
    def next(self):
        if self.itr >= len(self):
            raise StopIteration
        self.itr += 1
        return self.__getonerec(self.itr-1)
  
    def __iter__(self):
        self.itr = 0
        return self
  
    def intersect(self,chrom,querystart,querystop):
        """
        Fetch a intersection list with a given interval
        Return a list with all records within the interval
        dtype: numpy record array of (string,int,int,float)
        """
        intersectList = []
        if chrom in self.data:
            i = bisect.bisect(self.data[chrom]['end'],querystart)
            while (i < len(self.data[chrom])) and (self.data[chrom][i]['start'] < querystop):
                intersectList.append(self.data[chrom][i])
                i += 1
        if len(intersectList) == 0:
            return None
        else:
            intersectList = np.core.records.fromrecords(intersectList,dtype=self._bedgraphdtype)
            if intersectList[0]['start'] < querystart:
                intersectList[0]['start'] = querystart
            if intersectList[-1]['end'] > querystop:
                intersectList[-1]['end'] = querystop
            return intersectList
  
  #=========================================================
    def __getitem__(self,key):
        if isinstance(key,int):
            """For case a[1] or a[1:10], return a list of records"""
            return self.__getonerec(key)
        elif isinstance(key,slice):
            if key.step is None: step = 1
            else: step = key.step
            if key.start is None: start = 0
            else: start = key.start
            if key.stop is None: stop = len(self)
            else: stop = key.stop

            if start < 0: start += len(self)
            if stop < 0: stop += len(self)
            if start > len(self) or stop > len(self) :  raise IndexError, "The index out of range"
            records = []
            for i in range(start,stop,step):
                records.append(self.__getonerec(i))
            return np.core.records.fromrecords(records,dtype=self._bedgraphdtype)
      
        elif isinstance(key,tuple):
            """For case a['chr1',3000000:4000000], output average value"""
            chrom = key[0]
            if not chrom in self.data:
                raise KeyError, "Key %s doesn't exist!" % chrom
            try: 
                query = key[1]
            except Exception:
                raise TypeError, "Invalid argument type"
      
            assert isinstance(chrom,str)
            assert isinstance(query,slice)
            if query.start == None: 
                querystart = self.data[chrom][0]['start']
            else: 
                querystart = query.start
            if query.stop == None: 
                querystop = self.data[chrom][-1]['end']
            else: 
                querystop = query.stop
            assert querystop > querystart
            """
            Fetch all intersection with query.start:query.stop 
            """
            queryList = self.intersect(chrom,querystart,querystop)
            value = 0.0
            if not queryList is None:
                #Sum all value
                for rec in queryList:
                    value += rec['value'] * (rec['end'] - rec['start'])
    
            return value/(querystop-querystart)
        else:
            raise TypeError, "Invalid argument type"
    
    #=========================================================
    def __setitem__(self,key,value):
        assert isinstance(key,tuple)
        """For case a['chr1',3000000:4000000], input average value"""
        chrom = key[0]
        try: 
            query = key[1]
        except Exception:
            raise TypeError, "Invalid argument type"
        assert isinstance(chrom,str)
        assert isinstance(query,slice)
    
        new = np.array([(chrom,query.start,query.stop,value,'')],
                       dtype=self._bedgraphdtype)
        if not chrom in self.data:
            self.data[chrom] = []
            self.data[chrom].append(new)
            self.data[chrom] = np.core.records.fromrecords(self.data[chrom],dtype = self._bedgraphdtype)
            self._flush()
        else:
            i = bisect.bisect(self.data[chrom]['end'],query.start)
            deletelist = []
      
            if self.data[chrom][i]['start'] < query.start:
                self.data[chrom][i]['end'] = query.start
                i += 1
      
            insertLoc = i
            while (i < len(self.data[chrom])) and (self.data[chrom][i]['end'] < query.stop):
                #print self.data[chrom][i]
                deletelist.append(i)
                i += 1

            if i < len(self.data[chrom]):
                if self.data[chrom][i]['start'] < query.stop:
                    self.data[chrom][i]['start'] = query.stop
          
            self.data[chrom] = np.delete(self.data[chrom],deletelist)
            self.data[chrom] = np.insert(self.data[chrom],insertLoc,new)
        #add new data finished
    #=======================================================
    def filter(self,pattern):
        regpattern = re.compile(pattern)
        filterList = []
        for rec in self:
            if regpattern.match(rec['flag']):
                filterList.append(rec)
        return np.core.records.fromrecords(filterList,dtype=self._bedgraphdtype)
        
    #========================================================
    def save(self,filename,bedtype='bedgraph',style='%.8f'):
        """save bed file
        can be bedgraph,bedgraph with flag,bed
        """
        f = open(filename,'w')
        if bedtype == 'bedgraph':
            pattern = '%s\t%d\t%d\t'+style+'\n'
            for chrom in self.__sorted_keys:
                for line in self.data[chrom]:
                    f.write(pattern % (chrom,line['start'],line['end'],line['value']))
        elif bedtype == 'bed':
            pattern = '%s\t%d\t%d\t%s\n'
            for chrom in self.__sorted_keys:
                for line in self.data[chrom]:
                    f.write(pattern % (chrom,line['start'],line['end'],line['flag']))
        elif bedtype == 'bedgraph+':
            pattern = '%s\t%d\t%d\t'+style+'\t%s\n'
            for chrom in self.__sorted_keys:
                for line in self.data[chrom]:
                    f.write(pattern % (chrom,line['start'],line['end'],line['value'],line['flag']))
        elif bedtype == 'bed+':
            pattern = '%s\t%d\t%d\t%s\t'+style+'\n'
            for chrom in self.__sorted_keys:
                for line in self.data[chrom]:
                    f.write(pattern % (chrom,line['start'],line['end'],line['flag'],line['value']))
        else:
            raise TypeError, "Invalid argument type %s" % (bedtype)
        f.close

class modelgroup(object):
    def __init__(self,grouphandler,genome,idx):
        try:
            self.xyz = grouphandler['xyz'][:]
            self.r   = grouphandler['r'][:]
            self.log = cPickle.loads(grouphandler['log'].value)
            self.pym = cPickle.loads(grouphandler['pym'].value)
        except Exception as ex:
            raise ex
        self.score  = float(re.findall('Final score (\d+.\d+)',self.log)[0])
        self.genome = genome
        self.idx    = idx

        findRes = re.findall('#of Intra restraints: (\d+) #of Inter restraints: (\d+)',self.log)
        self.intraRestraints = float(findRes[0][0])
        self.interRestraints = float(findRes[0][1])
        findVio = re.findall('(\d+) violations in total',self.log)
        self.consecutiveViolations = float(findVio[0])
        try:
            self.contactViolations = float(findVio[1])            
        except:
            self.contactViolations = 0
            pass
    def __repr__(self):
        return 'Final Score: '+str(self.score)
    
    #-
    def getContactMap(self,contactRange=1):
        from scipy.spatial import distance
        n    = len(self.xyz)
        dist = distance.cdist(self.xyz,self.xyz,'euclidean')
        dcap = distance.cdist(self.r,-self.r)
        cmap = dist <= dcap*(contactRange+1)
        return cmap
    #=
    def getChromosomeCoordinates(self,chrom):
        Aids = np.flatnonzero(self.idx['chrom'] == chrom)
        Bids = Aids + len(self.idx)
        return self.xyz[Aids],self.xyz[Bids],self.r[Aids]
    #=
    
    def savepym(self,filename):
        pymfile = open(filename,'w')
        pymfile.write(self.pym)
        pymfile.flush()
        pymfile.close()
        return 0
    #=
    
#-
class modelstructures(object):
    """
    Instance manipulating hdf5 packed model structures
    take .hms file generated by tad modeling
    parameters:
    filename:  model result file *.hms
    usegrp:    group label array, usually assigned with probability prefix like ['p100','p020']
    """
    def __init__(self,filename,usegrp):
        if not os.path.isfile(filename):
            raise RuntimeError,"File %s doesn't exist!\n" % (filename)
        try:
            h5f = h5py.File(filename,'r')
        except:
            raise RuntimeError, "Invalid filetype, hdf5 file is required!"
        try:
            self.genome = cPickle.loads(h5f['genome'].value)
            self.idx    = h5f['idx'][:]
        except:
            raise RuntimeError, "Invalid file, genome or idx not found!"
        
        self._structs = []
        self.grpnames = usegrp
        for prefix in usegrp:
            try:
                grp = h5f[prefix]
            except:
                raise RuntimeError, "Group name %s doesn't exist!"%(prefix)
            mg = modelgroup(grp,self.genome,self.idx)
            self._structs.append(mg)
        #--
    #++++++++++++++++++++++++++++++++++++++++++++
    def __len__(self):
        return len(self._structs)
    
    def next(self):
        if self._itr >= len(self):
            raise StopIteration
        self._itr += 1
        return self._structs[self._itr-1]
  
    def __iter__(self):
        self._itr = 0
        return self
    
    def __getitem__(self,key):
        if isinstance(key,str):
            return self._structs[self.grpnames.index(key)]
        else:
            return self._structs[key]
    
        
if __name__=='__main__':
    pass
